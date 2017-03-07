/*
* Contar cuerpos celestes
*
* Asignatura Computación Paralela (Grado Ingeniería Informática)
* Código secuencial base
*
* @author Ana Moretón Fernández
* @version v1.2
*
* (c) 2017, Grupo Trasgo, Universidad de Valladolid
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "cputils.h"
#include <omp.h>


/* Substituir min por el operador */
#define min(x,y)    ((x) < (y)? (x) : (y))

/**
* Funcion secuencial para la busqueda de mi bloque
*/
int computation(int x, int y, int columns, int* matrixData, int *matrixResult, int *matrixResultCopy){
	// Inicialmente cojo mi indice
	int result=matrixResultCopy[x*columns+y];
	if( result!= -1){
		//Si es de mi mismo grupo, entonces actualizo
		if(matrixData[(x-1)*columns+y] == matrixData[x*columns+y])
		{
			result = min (result, matrixResultCopy[(x-1)*columns+y]);
		}
		if(matrixData[(x+1)*columns+y] == matrixData[x*columns+y])
		{
			result = min (result, matrixResultCopy[(x+1)*columns+y]);
		}
		if(matrixData[x*columns+y-1] == matrixData[x*columns+y])
		{
			result = min (result, matrixResultCopy[x*columns+y-1]);
		}
		if(matrixData[x*columns+y+1] == matrixData[x*columns+y])
		{
			result = min (result, matrixResultCopy[x*columns+y+1]);
		}

		// Si el indice no ha cambiado retorna 0
		if(matrixResult[x*columns+y] == result){ return 0; }
		// Si el indice cambia, actualizo matrix de resultados con el indice adecuado y retorno 1
		else { matrixResult[x*columns+y]=result; return 1;}

	}
	return 0;
}

/**
* Funcion principal
*/
int main (int argc, char* argv[])
{

	/* 1. Leer argumento y declaraciones */
	if (argc < 2) 	{
		printf("Uso: %s <imagen_a_procesar>\n", argv[0]);
		return(EXIT_SUCCESS);
	}
	char* image_filename = argv[1];

	int rows=-1;
	int columns =-1;
	int *matrixData=NULL;
	int *matrixResult=NULL;
	int *matrixResultCopy=NULL;
	int numBlocks=-1;



	/* 2. Leer Fichero de entrada e inicializar datos */

	/* 2.1 Abrir fichero */
	FILE *f = cp_abrir_fichero(image_filename);

	// Compruebo que no ha habido errores
	if (f==NULL)
	{
	   perror ("Error al abrir fichero.txt");
	   return -1;
	}

	/* 2.2 Leo valores del fichero */
	int i,j,valor;
	fscanf (f, "%d\n", &rows);
	fscanf (f, "%d\n", &columns);
	// Añado dos filas y dos columnas mas para los bordes
	rows=rows+2;
	columns = columns+2;

	/* 2.3 Reservo la memoria necesaria para la matriz de datos */
	matrixData= (int *)malloc( rows*(columns) * sizeof(int) );
	if ( (matrixData == NULL)   ) {
 		perror ("Error reservando memoria");
	   	return -1;
	}

	/* 2.4 Inicializo matrices */
	for(i=0;i< rows; i++){
		for(j=0;j< columns; j++){
			matrixData[i*(columns)+j]=-1;
		}
	}
	/* 2.5 Relleno bordes de la matriz */
	for(i=1;i<rows-1;i++){
		matrixData[i*(columns)+0]=0;
		matrixData[i*(columns)+columns-1]=0;
	}
	for(i=1;i<columns-1;i++){
		matrixData[0*(columns)+i]=0;
		matrixData[(rows-1)*(columns)+i]=0;
	}
	/* 2.6 Relleno la matriz con los datos del fichero */
	for(i=1;i<rows-1;i++){
		for(j=1;j<columns-1;j++){
			fscanf (f, "%d\n", &matrixData[i*(columns)+j]);
		}
	}
	fclose(f);

	#ifdef WRITE
		printf("Inicializacion \n");
		for(i=0;i<rows;i++){
			for(j=0;j<columns;j++){
				printf ("%d\t", matrixData[i*(columns)+j]);
			}
			printf("\n");
		}
	#endif


	/* PUNTO DE INICIO MEDIDA DE TIEMPO */
	double t_ini = cp_Wtime();

//
// EL CODIGO A PARALELIZAR COMIENZA AQUI
//

	/* 3. Etiquetado inicial */
	matrixResult= (int *)malloc( (rows)*(columns) * sizeof(int) );
	matrixResultCopy= (int *)malloc( (rows)*(columns) * sizeof(int) );
	if ( (matrixResult == NULL)  || (matrixResultCopy == NULL)  ) {
 		perror ("Error reservando memoria");
	   	return -1;
	}
	#pragma omp parallel for shared(matrixResult, matrixData) private(i,j)
	for(i=0;i< rows; i++){
		for(j=0;j< columns; j++){
			matrixResult[i*(columns)+j]=-1;
			// Si es 0 se trata del fondo y no lo computamos
			if(matrixData[i*(columns)+j]!=0){
				matrixResult[i*(columns)+j]=i*(columns)+j;
			}
		}
	}



	/* 4. Computacion */
	int t=0;
	/* 4.1 Flag para ver si ha habido cambios y si se continua la ejecucion */
	int flagCambio=1;

	// int sizeChunk = columns*rows/omp_get_max_threads();
	int sizeChunk = 1;
	int vRegions, hRegions, vBlocks, hBlocks, ri, rj;

	// Tiling size definition
	// - By number of regions
	vRegions = hRegions = omp_get_num_threads()*4;
	vBlocks = (rows-2)/vRegions, hBlocks = (columns-2)/hRegions;
	// - By size of block
	// vBlocks = hBlocks = omp_get_num_threads();
	// vRegions = (rows-2)/vBlocks, hRegions = (columns-2)/hBlocks;

	/* 4.2 Busqueda de los bloques similiares */
	for(t=0; flagCambio !=0; t++){
		flagCambio=0;

		/* 4.2.1 Actualizacion copia */
		// #pragma omp parallel shared(rows,columns,matrixData,matrixResult,matrixResultCopy,vBlocks,hBlocks,flagCambio)
		#pragma omp parallel shared(rows,columns,sizeChunk,vRegions,hRegions,vBlocks,hBlocks,flagCambio) firstprivate(matrixData,matrixResult,matrixResultCopy)
		{

			#pragma omp for private(ri,rj,i,j) schedule(guided,sizeChunk) collapse(2)
			for(ri=0;ri<=vRegions;ri++)
				for(rj=0;rj<=hRegions;rj++)
					for(i=(vBlocks*ri)+1;i<min(rows-1,(vBlocks*ri+vBlocks)+1);i++){
						for(j=(hBlocks*rj)+1;j<min(columns-1,(hBlocks*rj+hBlocks)+1);j++){
							if(matrixResult[i*(columns)+j]!=-1){
								matrixResultCopy[i*(columns)+j]=matrixResult[i*(columns)+j];
							}
						}
					}
			// #pragma omp for private(i,j) schedule(dynamic,sizeChunk)
			// for(i=1;i<rows-1;i++){
			// 	for(j=1;j<columns-1;j++){
			// 		if(matrixResult[i*(columns)+j]!=-1){
			// 			matrixResultCopy[i*(columns)+j]=matrixResult[i*(columns)+j];
			// 		}
			// 	}
			// }

			/* 4.2.2 Computo y detecto si ha habido cambios */
			#pragma omp for private(ri,rj,i,j) reduction(+:flagCambio) schedule(guided,sizeChunk) collapse(2)
			for(ri=0;ri<=vRegions;ri++)
				for(rj=0;rj<=hRegions;rj++)
					for(i=(vBlocks*ri)+1;i<min(rows-1,(vBlocks*ri+vBlocks)+1);i++){
						for(j=(hBlocks*rj)+1;j<min(columns-1,(hBlocks*rj+hBlocks)+1);j++){
							flagCambio= flagCambio+ computation(i,j,columns, matrixData, matrixResult, matrixResultCopy);
						}
					}
			// #pragma omp for private(i,j) reduction(+:flagCambio) schedule(dynamic,sizeChunk)
			// for(i=1;i<rows-1;i++){
			// 	for(j=1;j<columns-1;j++){
			// 		flagCambio= flagCambio+ computation(i,j,columns, matrixData, matrixResult, matrixResultCopy);
			// 	}
			// }
		}

		#ifdef DEBUG
			printf("\nResultados iter %d: \n", t);
			for(i=0;i<rows;i++){
				for(j=0;j<columns;j++){
					printf ("%d\t", matrixResult[i*columns+j]);
				}
				printf("\n");
			}
		#endif

	}

	/* 4.3 Inicio cuenta del numero de bloques */
	numBlocks=0;
	#pragma omp parallel for shared(matrixResult,rows,columns) private(i,j) reduction(+:numBlocks)
	for(i=1;i<rows-1;i++){
		for(j=1;j<columns-1;j++){
			if(matrixResult[i*columns+j] == i*columns+j) numBlocks++;
		}
	}

//
// EL CODIGO A PARALELIZAR TERMINA AQUI
//

	/* PUNTO DE FINAL DE MEDIDA DE TIEMPO */
 	double t_fin = cp_Wtime();



	/* 5. Comprobación de resultados */
  	double t_total = (double)(t_fin - t_ini);

	printf("Result: %d\n", numBlocks);
	printf("Time: %lf\n", t_total);
	#ifdef WRITE
		printf("Resultado: \n");
		for(i=0;i<rows;i++){
			for(j=0;j<columns;j++){
				printf ("%d\t", matrixResult[i*columns+j]);
			}
			printf("\n");
		}
	#endif

	/* 6. Liberacion de memoria */
	free(matrixData);
	free(matrixResult);
	free(matrixResultCopy);

}
