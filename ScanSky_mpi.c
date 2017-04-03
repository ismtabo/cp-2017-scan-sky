/*
* Contar cuerpos celestes
*
* Rodríguez Ares, Silvia
* Taboada Rodero, Ismael
*/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "cputils.h"
#include <mpi.h>


/* Substituir min por el operador */
#define min(x,y)    ((x) < (y)? (x) : (y))


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
	int world_rank = -1;
	int world_size = -1;
	double t_ini;
	int i,j;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &world_size);

	if ( world_rank == 0 ) {

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
		int valor;
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


		/* PUNTO DE INICIO MEDIDA DE TIEMPO */
		t_ini = cp_Wtime();
	}

	//
	// EL CODIGO A PARALELIZAR COMIENZA AQUI
	//

    int k;
    int dimensions[3];  // Vector of matrix and submatrix dimensions for communication
    int size_block, ncells;  // Number of rows and number of cells of proc submatrix
    int displacement; // Displacement of proc at initial matrix
    int vectorSizeBlocks[world_size], vectorDis[world_size];  // Vector with number of cells of and displacement for each proc
    int previous, next;  // Rank of previos and next proc 
    // TODO: Add submatrix structure for matrixData
    int *sub_matrixResultCopy, *sub_matrixResult;  // Proper submatrixes 
    MPI_Status stat;
    MPI_Request req;

	if ( world_rank == 0 ) {

		/* 3. Etiquetado inicial */
		matrixResult= (int *)malloc( (rows)*(columns) * sizeof(int) );
		matrixResultCopy= (int *)malloc( (rows)*(columns) * sizeof(int) );
		if ( (matrixResult == NULL)  || (matrixResultCopy == NULL)  ) {
			perror ("Error reservando memoria");
			return -1;
		}
		for(i=0;i< rows; i++){
			for(j=0;j< columns; j++){
				matrixResult[i*(columns)+j]=-1;
				// Si es 0 se trata del fondo y no lo computamos
				if(matrixData[i*(columns)+j]!=0){
					matrixResult[i*(columns)+j]=i*(columns)+j;
				}
			}
		}
        dimensions[0] = rows;
        dimensions[1] = columns;
        // TODO: Change function to balance deal
        dimensions[2] = size_block = rows/world_size; 

    }
    #ifdef DEBUG
    if(world_rank==0){
        printf("Matrix\n");
        for(i=0; i<rows; ++i){
            for(j=0; j<columns; ++j){
                printf("%d ", matrixResult[i*columns+j]);
            }
            printf("\n");
        }
        fflush(stdout);
    }
    #endif

    // Broadcast matrix dimensions
    MPI_Bcast(dimensions, 3, MPI_INT, 0, MPI_COMM_WORLD);
    previous = world_rank - 1; 
    next = world_rank + 1; 

    if( world_rank ){
        rows = dimensions[0];
        columns = dimensions[1];
        size_block = dimensions[2];
        matrixData = (int*) malloc(sizeof(int)*rows*columns);
    }
    
    if( world_rank == world_size-1)
        size_block += rows - size_block*world_size;
    
    // Broadcast of data matrix
    // TODO: Change Bcast into scatter of:
    // - proper submatrix
    // - first rows for each process
    // - last row for each process
    MPI_Bcast(matrixData, rows*columns, MPI_INT, 0, MPI_COMM_WORLD);

    ncells = size_block*columns;

    // Group number of cells for each proc submatrix
    MPI_Gather(&ncells, 1, MPI_INT, vectorSizeBlocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate displacements
    if(world_rank==0){
        vectorDis[0] = 0;
        for(i=1; i<world_size; i++)
            // TODO: Calculate displacements for:
            // - Scatter of first rows
            // - Scatter of last rows
            vectorDis[i] = vectorDis[i-1]+vectorSizeBlocks[i-1];

        #ifdef DEBUG
        for(i=0; i<world_size; i++)
            printf("[%d] %d %d ",i, vectorSizeBlocks[i], vectorDis[i]);
        printf("\n");fflush(stdout);
        #endif
    }

    // Scatter of each displacement
    MPI_Scatter(vectorDis, 1, MPI_INT, &displacement, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(world_size > 1){
        if( world_rank == 0 ){
            size_block +=1;
            previous = MPI_PROC_NULL;
        } else if( world_rank < world_size-1 ) {
            size_block += 2;
        } else if( world_rank == world_size-1) {
            size_block += 1;
            next = MPI_PROC_NULL; 
        }
    }

    sub_matrixResult = (int*) malloc(sizeof(int)*size_block*columns);
    sub_matrixResultCopy = (int*) malloc(sizeof(int)*size_block*columns);
    for(i=0; i<size_block; ++i)
        for(j=0; j<columns; ++j){
            sub_matrixResult[i*columns+j]=-1;
            sub_matrixResultCopy[i*columns+j]=-1;
        }
            
    
    MPI_Scatterv(matrixResult, vectorSizeBlocks, vectorDis, MPI_INT, sub_matrixResult+(world_rank?columns:0), ncells, MPI_INT, 0, MPI_COMM_WORLD);
    
    #ifdef DEBUG
    if(world_rank==4){
         printf("[%d] Init submatrix\n", world_rank);
         for(i=0; i<size_block; ++i){
             for(j=0; j<columns; ++j){
                 printf("%d ", sub_matrixResult[i*columns+j]);
             }
             printf("\n");
         }
        fflush(stdout);
    }
    #endif

	/* 4. Computacion */
	int t=0;
	/* 4.1 Flag para ver si ha habido cambios y si se continua la ejecucion */
	int flagCambio=1; 

	/* 4.2 Busqueda de los bloques similiares */
	for(t=0; flagCambio !=0; t++){
		flagCambio=0; 
    
        if(world_size > 1){
            // Segment fault at sending on both programs
            MPI_Isend(&sub_matrixResult[columns], columns, MPI_INT, previous, 0, MPI_COMM_WORLD, &req);  // Send row #1
            MPI_Recv(&sub_matrixResult[(size_block-1)*columns], columns, MPI_INT, next, 0, MPI_COMM_WORLD, &stat); // Wait to last row
            MPI_Isend(&sub_matrixResult[(size_block-2)*columns], columns, MPI_INT, next, 0, MPI_COMM_WORLD, &req);  // Send penultimate row1
            MPI_Recv(&sub_matrixResult[0], columns, MPI_INT, previous, 0, MPI_COMM_WORLD, &stat); // Wait to first row

            #ifdef DEBUG
            if(world_rank==0){
                printf("[%d]%d Submatrix\n",world_rank, t);
                for(i=0; i<size_block; ++i){
                    for(j=0; j<columns; ++j){
                        printf("%d ", sub_matrixResult[i*columns+j]);
                    }
                    printf("\n");
                }
                fflush(stdout);
            }
            #endif
        }

		/* 4.2.1 Actualizacion copia */
		for(i=0;i<size_block;i++){
			for(j=0;j<columns;j++){
				if(sub_matrixResult[i*(columns)+j]!=-1){
					sub_matrixResultCopy[i*(columns)+j]=sub_matrixResult[i*(columns)+j];
				}
			}
		}
        #ifdef DEBUG
        if(world_rank==0){
            printf("[%d]%d Submatrix copy\n",world_rank, t);
            for(i=0; i<size_block; ++i){
                for(j=0; j<columns; ++j){
                    printf("%d ", sub_matrixResultCopy[i*columns+j]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        #endif

		/* 4.2.2 Computo y detecto si ha habido cambios */
		for(i=1;i<size_block-1;i++){
			for(k=1;k<columns-1;k++){
    			int result,sol;
                j=i*columns+k;
                int dataj = (world_rank==0?i:i-1)*columns+k+displacement;
                result = sub_matrixResultCopy[j];
    			sol=0;
                if(result!=-1){
    		    	//Si es de mi mismo grupo, entonces actualizo
    		        if(matrixData[dataj-columns] == matrixData[dataj])
    		        {
    		            result = min (result, sub_matrixResultCopy[j-columns]);
    		        }
    		        if(matrixData[dataj+columns] == matrixData[dataj])
    		        {
    		            result = min (result, sub_matrixResultCopy[j+columns]);
    		        }
    		        if(matrixData[dataj-1] == matrixData[dataj])
    		        {
    		            result = min (result, sub_matrixResultCopy[j-1]);
    		        }
    		        if(matrixData[dataj+1] == matrixData[dataj])
    		        {
    		            result = min (result, sub_matrixResultCopy[j+1]);
    		        }
    		        // Si el indice no ha cambiado retorna 0
    		        if(sub_matrixResult[j] == result){ sol=0; }
    		        // Si el indice cambia, actualizo sub_matrix de resultados con el indice adecuado y retorno 1
    		        else { sub_matrixResult[j]=result; sol=1;}
    		        flagCambio= flagCambio+ sol;
                }
		    }
		}

        #ifdef DEBUG
        if(world_rank==4){
            printf("[%d]%d Submatrix result\n",world_rank, t);
            for(i=0; i<size_block; ++i){
                for(j=0; j<columns; ++j){
                    printf("%d ", sub_matrixResult[i*columns+j]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
        #endif

        MPI_Allreduce(MPI_IN_PLACE, &flagCambio, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}

    // MPI_Gatherv(sub_matrixResult+(world_rank?columns:0), ncells, MPI_INT, matrixResult, vectorSizeBlocks, vectorDis, MPI_INT, 0, MPI_COMM_WORLD);   

    // if(world_rank==0)
    // {
	// 	/* 4.3 Inicio cuenta del numero de bloques */
    //     printf("Matrix result\n");
    // 	for(i=0;i<rows;i++){
    // 		for(j=0;j<columns;j++){
    // 			printf ("%d\t", matrixResult[i*columns+j]);
    // 		}
    // 		printf("\n");
    // 	}
	// 	numBlocks=0;
	// 	for(i=1;i<rows-1;i++){
	// 		for(j=1;j<columns-1;j++){
	// 			if(matrixResult[i*columns+j] == i*columns+j) numBlocks++; 
	// 		}
	// 	}
	// }

    // Count blobs on subMatrix
    int nBlocks=0;
    for(i=1;i<size_block-1;i++){
    	for(j=1;j<columns-1;j++){
    		if(sub_matrixResult[i*columns+j] == (world_rank==0?i:i-1)*columns+j+displacement) nBlocks++; 
    	}
    }        
    
    MPI_Reduce(&nBlocks, &numBlocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(world_rank)
        free(matrixData);    
    free(sub_matrixResult);
    free(sub_matrixResultCopy);

	//
	// EL CODIGO A PARALELIZAR TERMINA AQUI
	//
	if ( world_rank == 0 ) {

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

	MPI_Finalize();
	return 0;
}
