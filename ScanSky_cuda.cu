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
#define T 1024;
/**
* Funcion paralela para inicializar la matriz de etiquetas
*/
__global__ void init_matrixResult(int rows, int columns, int *matrixData, int *matrixResult, int *matrixResultCopy, int *flagCambioArray) {
	int i = blockDim.y * blockIdx.y + threadIdx.y;
	int j = blockDim.x * blockIdx.x + threadIdx.x;

	if (0 <= i && i <= rows - 1 && 0 <= j && j <= columns - 1){
		matrixResultCopy[i*(columns)+j]=-1;
		matrixResult[i*(columns)+j]=-1;
		// Si es 0 se trata del fondo y no lo computamos
		if(matrixData[i*(columns)+j]!=0){
			matrixResult[i*(columns)+j]=i*(columns)+j;
		}
		flagCambioArray[i*(columns)+j]=0;
	}
}

/**
* Funcion paralela para la copia de matrixResult a matrixResultCopy
*/
__global__ void cpy_matrixResult(int rows, int columns, int* matrixResult, int* matrixResultCopy) {
	int i = blockDim.y * blockIdx.y + threadIdx.y;
	int j = blockDim.x * blockIdx.x + threadIdx.x;

	if (0 <= i && i <= rows - 1 && 0 <= j && j <= columns - 1){
		if(matrixResult[i*(columns)+j]!=-1){
			matrixResultCopy[i*(columns)+j]=matrixResult[i*(columns)+j];
		}
	}
}

/**
* Funcion paralela para la busqueda de mi bloque
*/
__global__ void computation(int rows, int columns, int* matrixData, int *matrixResult, int *matrixResultCopy, int *flagCambioArray){
	int x = blockDim.y * blockIdx.y + threadIdx.y;
	int y = blockDim.x * blockIdx.x + threadIdx.x;


	if (0 < x && x < rows - 1 && 0 < y && y < columns - 1){
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
			if(matrixResult[x*columns+y] == result){ flagCambioArray[x*columns+y] = 0; }
			// Si el indice cambia, actualizo matrix de resultados con el indice adecuado y retorno 1
			else { matrixResult[x*columns+y]=result; flagCambioArray[x*columns+y] = 1;}
		}
	}
}

/**
* Funcion secuencial para la busqueda de mi bloque
*/
int _computation(int x, int y, int columns, int* matrixData, int *matrixResult, int *matrixResultCopy){
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
//funcion paralela reduccion
__global__ void reduce(int rows, int columns,int *g_idata, int *g_odata) {

    extern __shared__ int sdata[];

    // each thread loads one element from global to shared mem
    // note use of 1D thread indices (only) in this kernel
		int i = blockIdx.y*blockDim.y + threadIdx.y;
    int j = blockIdx.x*blockDim.x + threadIdx.x;

		if(0 < i && i < rows-1 && 0 < j && j < columns-1){
			sdata[threadIdx.y*blockDim.x+threadIdx.x] = g_idata[i*columns+j];
		} else {
			sdata[threadIdx.y*blockDim.x+threadIdx.x] = 0;
		}

		__syncthreads();
		// do reduction in shared mem
		for (int s=1; s < blockDim.x*blockDim.y; s *=2)
		{
			int index = 2 * s * (threadIdx.y*blockDim.x + threadIdx.x);

			if (index < blockDim.x*blockDim.y)
			{
				sdata[index] += sdata[index + s];
			}
			__syncthreads();
		}

		// write result for this block to global mem
		if (threadIdx.y==0 && threadIdx.x == 0)
			atomicAdd(g_odata,sdata[0]);

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
	cudaError_t error;

	// Vector auxiliar para la persistencia de flagCambio de cada hilo
	int *flagCambioA = NULL;
	int *flagCambioB = NULL;

	// Vectores auxiliares utilizados en dispositivo
	int *matrixDataDev = NULL;
	int *matrixResultDev = NULL;
	int *matrixResultCopyDev = NULL;
	int *flagCambioADev = NULL;
	int *flagCambioBDev = NULL;

	cudaMalloc((void **) &matrixDataDev, (rows) * (columns) * sizeof(int));
	cudaMalloc((void **) &matrixResultDev, (rows) * (columns) * sizeof(int));
	cudaMalloc((void **) &matrixResultCopyDev, (rows) * (columns) * sizeof(int));
	cudaMalloc((void **) &flagCambioADev, (rows) * (columns) * sizeof(int));
	cudaMalloc((void **) &flagCambioBDev, (rows) * (columns) * sizeof(int));
	/* 3. Etiquetado inicial */
	matrixResult= (int *)malloc( (rows)*(columns) * sizeof(int) );
	matrixResultCopy= (int *)malloc( (rows)*(columns) * sizeof(int) );
	flagCambioA = (int *)malloc( (rows)*(columns) * sizeof(int) );
	flagCambioB = (int *)malloc( (rows)*(columns) * sizeof(int) );
	flagCambioB[0]=0;
	if ( (matrixResult == NULL)  || (matrixResultCopy == NULL) || (flagCambioA == NULL) || (flagCambioB == NULL))  {
 		perror ("Error reservando memoria");
	   	return -1;
	}
// Inicializacion de matrices paralelizada por init_matrixResult
	// for(i=0;i< rows; i++){
	// 	for(j=0;j< columns; j++){
	// 		matrixResultCopy[i*(columns)+j]=-1;
	// 		matrixResult[i*(columns)+j]=-1;
	// 		// Si es 0 se trata del fondo y no lo computamos
	// 		if(matrixData[i*(columns)+j]!=0){
	// 			matrixResult[i*(columns)+j]=i*(columns)+j;
	// 		}
	// 	}
	// }
	int blockDim_x, blockDim_y, gridDim_x, gridDim_y;
	blockDim_x = 32; blockDim_y = 8;
	gridDim_y = rows / blockDim_y + (rows % blockDim_y ? 1 : 0);
	gridDim_x = columns / blockDim_x + (columns % blockDim_x ? 1 : 0);
	dim3 blockDims(blockDim_x, blockDim_y);
	dim3 gridDims(gridDim_x, gridDim_y);

	error = cudaMemcpy(matrixDataDev, matrixData, (rows) * (columns) * sizeof(int), cudaMemcpyHostToDevice);
	if ( error != cudaSuccess )
		printf("ErrCUDA cpy mD h2d: %s\n", cudaGetErrorString( error ) );

	// Lanzamiento del kernel
	init_matrixResult<<<gridDims, blockDims>>>(rows, columns, matrixDataDev, matrixResultDev, matrixResultCopyDev, flagCambioADev);
	error = cudaGetLastError();
	if ( error != cudaSuccess )
		printf("ErrCUDA init mR d: %s\n", cudaGetErrorString( error ) );

	// error = cudaMemcpy(matrixResult, matrixResultDev, (rows) * (columns) * sizeof(int), cudaMemcpyDeviceToHost);
	// if ( error != cudaSuccess )
	// 	printf("ErrCUDA cpy mR d2h: %s\n", cudaGetErrorString( error ) );
	// error = cudaMemcpy(matrixResultCopy, matrixResultCopyDev, (rows) * (columns) * sizeof(int), cudaMemcpyDeviceToHost);
	// if ( error != cudaSuccess )
	// 	printf("ErrCUDA cpy mRC d2h: %s\n", cudaGetErrorString( error ) );

	/* 4. Computacion */
	int t=0;
	/* 4.1 Flag para ver si ha habido cambios y si se continua la ejecucion */
	int flagCambio=1;

	/* 4.2 Busqueda de los bloques similiares */
	for(t=0; flagCambio !=0; t++){
		flagCambio=0;
		flagCambioB[0]=0;
		cudaMemcpy(flagCambioBDev, flagCambioB, sizeof(int), cudaMemcpyHostToDevice);

		/* 4.2.1 Actualizacion copia */
// Copia de matrices paralelizada por cpy_matrixResult
		// for(i=1;i<rows-1;i++){
		// 	for(j=1;j<columns-1;j++){
		// 		if(matrixResult[i*(columns)+j]!=-1){
		// 			matrixResultCopy[i*(columns)+j]=matrixResult[i*(columns)+j];
		// 		}
		// 	}
		// }
		// error = cudaMemcpy(matrixResultDev, matrixResult, (rows) * (columns) * sizeof(int), cudaMemcpyHostToDevice);
		// if ( error != cudaSuccess )
		// 	printf("ErrCUDA cpy mR d2h: %s\n", cudaGetErrorString( error ) );

		cpy_matrixResult<<<gridDims, blockDims>>>(rows, columns, matrixResultDev, matrixResultCopyDev);
		error = cudaGetLastError();
		if ( error != cudaSuccess )
			printf("ErrCUDA cpy mR d: %s\n", cudaGetErrorString( error ) );

		// error = cudaMemcpy(matrixResultCopy, matrixResultCopyDev, (rows) * (columns) * sizeof(int), cudaMemcpyDeviceToHost);
		// if ( error != cudaSuccess )
		// 	printf("ErrCUDA cpy mRC d2h: %s\n", cudaGetErrorString( error ) );

		/* 4.2.2 Computo y detecto si ha habido cambios */
// Computo paralelizado por computation
		// for(i=1;i<rows-1;i++){
		// 	for(j=1;j<columns-1;j++){
		// 		flagCambio= flagCambio+ computation(i,j,columns, matrixData, matrixResult, matrixResultCopy);
		// 	}
		// }
		computation<<<gridDims, blockDims>>>(rows, columns, matrixDataDev, matrixResultDev, matrixResultCopyDev, flagCambioADev);
		error = cudaGetLastError();
		if ( error != cudaSuccess )
			printf("ErrCUDA computation d: %s\n", cudaGetErrorString( error ) );


		reduce<<<gridDims, blockDims, 256*sizeof(int)>>>(rows, columns, flagCambioADev, flagCambioBDev);
		error = cudaGetLastError();
		if ( error != cudaSuccess )
			printf("ErrCUDA reduce d: %s\n", cudaGetErrorString( error ) );


		error = cudaMemcpy(flagCambioB, flagCambioBDev, sizeof(int), cudaMemcpyDeviceToHost);
		if ( error != cudaSuccess )
			printf("ErrCUDA cpy fCA d2h: %s\n", cudaGetErrorString( error ) );

/*
		error = cudaMemcpy(flagCambioA, flagCambioADev, rows*columns*sizeof(int), cudaMemcpyDeviceToHost);
		if ( error != cudaSuccess )
			printf("ErrCUDA cpy fCA d2h: %s\n", cudaGetErrorString( error ) );

		for(i=0; i<(rows)*(columns); i++){
			flagCambio = flagCambio + flagCambioA[i];
		}
*/
		 flagCambio=flagCambioB[0];
		//flagCambio=0;
		// printf("flagCambio: %d\n", flagCambioB[0]);

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

	// Recuperación de matrixResult en anfitrion
	error = cudaMemcpy(matrixResult, matrixResultDev, (rows) * (columns) * sizeof(int), cudaMemcpyDeviceToHost);
	if ( error != cudaSuccess )
		printf("ErrCUDA fetch mR d2h: %s\n", cudaGetErrorString( error ) );

	/* 4.3 Inicio cuenta del numero de bloques */
	numBlocks=0;
	for(i=1;i<rows-1;i++){
		for(j=1;j<columns-1;j++){
			if(matrixResult[i*columns+j] == i*columns+j) numBlocks++;
		}
	}

	free(flagCambioA);
	free(flagCambioB);

	//cudaFree();
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
