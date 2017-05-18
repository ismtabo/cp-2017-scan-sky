/*
* Contar cuerpos celestes
*
* Asignatura Computación Paralela (Grado Ingeniería Informática)
* Código secuencial base
*
* @author Ana Moretón Fernández, Arturo Gonzalez-Escribano
* @version v1.3
*
* (c) 2017, Grupo Trasgo, Universidad de Valladolid
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <cuda.h>
#include "cputils.h"


/* Substituir min por el operador */
#define min(x,y)    ((x) < (y)? (x) : (y))
#define TSM 2048 // max threads per SM
#define T 1024 //max threads per block
#define BSM 16 //max blocks per SM
/**
* Funcion secuencial para la busqueda de mi bloque
*/
__global__ void update_copy(int rows, int columns, int *matrixResultDev, int *matrixResultCopyDev){
  int col, row;
	col = blockIdx.x*blockDim.x+threadIdx.x;
	row = blockIdx.y*blockDim.y+threadIdx.y;
	if(0 < row && row < rows && 0 < col && col < columns){
		if(matrixResultDev[row*(columns)+col]!=-1){
			matrixResultCopyDev[row*(columns)+col]=matrixResultDev[row*(columns)+col];
		}
	}
}

__global__ void computation(int rows, int columns, int* matrixData, int *matrixResult, int *matrixResultCopy, int *flagCambioArrayDev){
	//flagCambio= flagCambio+ computation(i,j,columns, matrixData, matrixResult, matrixResultCopy);
	//		computation<<<dimGrid, dimBlock>>>(columns,matrixDataDev,matrixResultDev,matrixResultCopyDev,flagCambioArrayDev);
	int col, row, flagCambio=0;
	col = blockIdx.x*blockDim.x+threadIdx.x;
	row = blockIdx.y*blockDim.y+threadIdx.y;
	// Inicialmente cojo mi indice
	if (0 < row && row < rows && 0 < col && col < columns) {
		int result=matrixResultCopy[row*columns+col];
		if( result!= -1){
			//Si es de mi mismo grupo, entonces actualizo
			if(matrixData[(row-1)*columns+col] == matrixData[row*columns+col])
			{
				result = min (result, matrixResultCopy[(row-1)*columns+col]);
			}
			if(matrixData[(row+1)*columns+col] == matrixData[row*columns+col])
			{
				result = min (result, matrixResultCopy[(row+1)*columns+col]);
			}
			if(matrixData[row*columns+col-1] == matrixData[row*columns+col])
			{
				result = min (result, matrixResultCopy[row*columns+col-1]);
			}
			if(matrixData[row*columns+col+1] == matrixData[row*columns+col])
			{
				result = min (result, matrixResultCopy[row*columns+col+1]);
			}

			// Si el indice no ha cambiado retorna 0
			if(matrixResult[row*columns+col] == result){ flagCambio=0; }
			// Si el indice cambia, actualizo matrix de resultados con el indice adecuado y retorno 1
			else { matrixResult[row*columns+col]=result; flagCambio=1;}

		}
	}
	flagCambioArrayDev[row*columns+col]=flagCambio;
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
	int i,j;
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

	cudaSetDevice(0);
	cudaDeviceSynchronize();

	/* PUNTO DE INICIO MEDIDA DE TIEMPO */
	double t_ini = cp_Wtime();

//
// EL CODIGO A PARALELIZAR COMIENZA AQUI
//

	/* 3. Etiquetado inicial */
  cudaError_t error;
	int *flagCambioArray;
	matrixResult= (int *)malloc( (rows)*(columns) * sizeof(int) );
	matrixResultCopy= (int *)malloc( (rows)*(columns) * sizeof(int) );
	flagCambioArray= (int *)malloc( (rows)*(columns) * sizeof(int) );

	if ( (matrixResult == NULL)  || (matrixResultCopy == NULL) || (flagCambioArray == NULL) ) {
 		perror ("Error reservando memoria");
	   	return -1;
	}
	for(i=0;i< rows; i++){
		for(j=0;j< columns; j++){
			matrixResultCopy[i*(columns)+j]=-1;
			matrixResult[i*(columns)+j]=-1;
			// Si es 0 se trata del fondo y no lo computamos
			if(matrixData[i*(columns)+j]!=0){
				matrixResult[i*(columns)+j]=i*(columns)+j;
			}
		}
	}

	int *matrixResultDev;
	int *matrixResultCopyDev;
	int *matrixDataDev;
	int *flagCambioArrayDev;

  // Reserva de memoria en dispositivo
  #ifdef ERROR
  	error = cudaMalloc((void **) &matrixResultDev, rows*(columns) * sizeof(int)  );
    if ( error != cudaSuccess )
      printf("ErrCUDA malloc mR: %s\n", cudaGetErrorString( error ) );
  	error = cudaMalloc((void **) &matrixResultCopyDev, rows*(columns) * sizeof(int) );
    if ( error != cudaSuccess )
      printf("ErrCUDA malloc mRC: %s\n", cudaGetErrorString( error ) );
  	error = cudaMalloc((void **) &matrixDataDev, rows*(columns) * sizeof(int)  );
    if ( error != cudaSuccess )
      printf("ErrCUDA malloc mD: %s\n", cudaGetErrorString( error ) );
  	error = cudaMalloc((void **) &flagCambioArrayDev, rows*(columns) * sizeof(int)  );
    if ( error != cudaSuccess )
      printf("ErrCUDA malloc fC: %s\n", cudaGetErrorString( error ) );
  #else
    cudaMalloc((void **) &matrixResultDev, rows*(columns) * sizeof(int)  );
    cudaMalloc((void **) &matrixResultCopyDev, rows*(columns) * sizeof(int) );
    cudaMalloc((void **) &matrixDataDev, rows*(columns) * sizeof(int)  );
    cudaMalloc((void **) &flagCambioArrayDev, rows*(columns) * sizeof(int)  );
  #endif

  // Copia de matrices de anfitrion a dispositivo
  #ifdef ERROR
  	error = cudaMemcpy(matrixResultDev,matrixResult,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
    if ( error != cudaSuccess )
      printf("ErrCUDA cpy mR: %s\n", cudaGetErrorString( error ) );
  	error = cudaMemcpy(matrixResultCopyDev,matrixResultCopy,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
    if ( error != cudaSuccess )
      printf("ErrCUDA cpy mRC: %s\n", cudaGetErrorString( error ) );
  	error = cudaMemcpy(matrixDataDev,matrixData,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
    if ( error != cudaSuccess )
      printf("ErrCUDA cpy mD: %s\n", cudaGetErrorString( error ) );
  #else
    cudaMemcpy(matrixResultDev,matrixResult,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(matrixResultCopyDev,matrixResultCopy,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(matrixDataDev,matrixData,rows*(columns) * sizeof(int), cudaMemcpyHostToDevice);
  #endif

	int rowblock, colblock;
	rowblock = colblock = 16;
	int rowgrid, colgrid;
	rowgrid = rows / rowblock + (rows % rowblock ? 1 : 0);
	colgrid = columns / colblock + (columns % colblock ? 1 : 0);
	dim3 dimGrid(rowgrid, colgrid);
	dim3 dimBlock(rowblock, colblock);


	/* 4. Computacion */
	int t=0;
	/* 4.1 Flag para ver si ha habido cambios y si se continua la ejecucion */
	int flagCambio=1;

	/* 4.2 Busqueda de los bloques similiares */
	for(t=0; flagCambio !=0; t++){
		flagCambio=0;

		/* 4.2.1 Actualizacion copia */
		/*
		for(i=1;i<rows-1;i++){
			for(j=1;j<columns-1;j++){
				if(matrixResult[i*(columns)+j]!=-1){
					matrixResultCopy[i*(columns)+j]=matrixResult[i*(columns)+j];
				}
			}
		}*/
		//kernel actualizacion copia
		update_copy<<<dimGrid, dimBlock>>>(rows, columns, matrixResultDev,matrixResultCopyDev);
    #ifdef ERROR
      error = cudaGetLastError();
      if ( error != cudaSuccess )
        printf("ErrCUDA update_copy: %s\n", cudaGetErrorString( error ) );
    #endif
		/* 4.2.2 Computo y detecto si ha habido cambios */
		/*
		for(i=1;i<rows-1;i++){
			for(j=1;j<columns-1;j++){
				flagCambio= flagCambio+ computation(i,j,columns, matrixData, matrixResult, matrixResultCopy);
			}
		}
		*/
		computation<<<dimGrid, dimBlock>>>(rows, columns,matrixDataDev,matrixResultDev,matrixResultCopyDev,flagCambioArrayDev);
    #ifdef ERROR
      error = cudaGetLastError();
      if ( error != cudaSuccess )
        printf("ErrCUDA computation: %s\n", cudaGetErrorString( error ) );
    #endif

    // Copia del array de flagCambio del dispositivo a anfitrion
    #ifdef ERROR
      error = cudaMemcpy(flagCambioArray,flagCambioArrayDev,rows*(columns) * sizeof(int), cudaMemcpyDeviceToHost);
      if ( error != cudaSuccess )
        printf("ErrCUDA cpy flagCambio: %s\n", cudaGetErrorString( error ) );
    #else
      cudaMemcpy(flagCambioArray,flagCambioArrayDev,rows*(columns) * sizeof(int), cudaMemcpyDeviceToHost);
    #endif

		for(i=0 ;i<rows*columns;i++){
			flagCambio=flagCambio+flagCambioArray[i];
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

  // Copia de matriz resultado de dispositivo a anfitrion
	cudaMemcpy(matrixResult,matrixResultDev,rows*(columns) * sizeof(int), cudaMemcpyDeviceToHost);

	/* 4.3 Inicio cuenta del numero de bloques */
	numBlocks=0;
	for(i=1;i<rows-1;i++){
		for(j=1;j<columns-1;j++){
			if(matrixResult[i*columns+j] == i*columns+j) numBlocks++;
		}
	}

  // Liberacion de estructuras creadas en dispositivo
	cudaFree(matrixResultDev);
	cudaFree(matrixResultCopyDev);
	cudaFree(matrixDataDev);
	cudaFree(flagCambioArrayDev);
  free(flagCambioArray);

//
// EL CODIGO A PARALELIZAR TERMINA AQUI
//

	/* PUNTO DE FINAL DE MEDIDA DE TIEMPO */
	cudaDeviceSynchronize();
 	double t_fin = cp_Wtime();


	/* 5. Comprobación de resultados */
  	double t_total = (double)(t_fin - t_ini);

	printf("Result: %d:%d\n", numBlocks, t);
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

	return 0;
}
