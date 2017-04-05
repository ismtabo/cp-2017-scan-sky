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
#include <string.h>


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
    int dimensions[2];  // Vector of matrix dimensions for communication
    int sub_rows, displacement, ncells;  // Number of rows, displacement and number of cells of proc
    int vectorRows[world_size], vectorSizes[world_size], vectorDis[world_size];  // Vector with number of rows, number of cells of and displacement for each proc
    int namelen, namelens[world_size];
    char name[MPI_MAX_PROCESSOR_NAME], names[world_size][MPI_MAX_PROCESSOR_NAME];
    int contIndex;
    int new_ranks[world_size];
    int previous, next;  // Rank of previos and next proc 
    int *sub_matrixData, *sub_matrixResult, *sub_matrixResultCopy;  // Proper submatrixes 
    MPI_Status stats[2];
    MPI_Request reqs[2];

    // Calculate previous and next procs' ranks
    previous = world_rank - 1; 
    next = world_rank + 1; 

    MPI_Get_processor_name(name, &namelen);
    MPI_Gather(&namelen, 1, MPI_INT, namelens, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(name, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, names, MPI_MAX_PROCESSOR_NAME, MPI_CHAR, 0, MPI_COMM_WORLD);

    if ( world_rank == 0) {
        // Group ranks by name of processor
        int start, next_start;
        char actual_proccesor[MPI_MAX_PROCESSOR_NAME];
        contIndex = 0;
        next_start = 0;
        strcpy(names[1], "ultron0x01");
        do {
            start = next_start;
            next_start = 0;
            new_ranks[contIndex++] = start;
            strcpy(actual_proccesor, names[start]);
#ifdef DEBUG
            printf("Actual: %s %d looking into %d\n", actual_proccesor, start, world_size);fflush(stdout);
#endif
            for(i=start+1; i<world_size; i++){
#ifdef DEBUG
                printf("CMP for %s %d\n", actual_proccesor, i);fflush(stdout);
#endif
                if(strcmp(names[i],"")!=0 && strcmp(actual_proccesor, names[i])==0){
                    new_ranks[contIndex++] = i;                     
#ifdef DEBUG
                    printf("Same %s %s\n", actual_proccesor, names[i]);fflush(stdout);
#endif
                    strcpy(names[i],""); 
                } else if(strcmp(names[i], "")!=0){
                    if(!next_start)
                        next_start = i;
#ifdef DEBUG
                    printf("Diff %s %s\n", actual_proccesor, names[i]);fflush(stdout);
#endif
                } 
#ifdef DEBUG
                printf("Compared %d\n", i);fflush(stdout);
#endif
            } 
#ifdef DEBUG
            printf("Found all for %s\n", actual_proccesor);fflush(stdout);
#endif
        } while(next_start);
//#ifdef DEBUG
        printf("New ranks: \n");
        for(i = 0; i < world_size; i++)
            printf("[%d]%d ", i, new_ranks[i]);
        printf("\n");fflush(stdout);
//#endif
    } 

    MPI_Bcast(new_ranks, world_size, MPI_INT, 0, MPI_COMM_WORLD);
    for(i = 0; world_rank!=new_ranks[i] && i < world_size; i++){}
    
    previous = (i-1 >= 0)? new_ranks[i-1] : MPI_PROC_NULL;
    next = (i+1 < world_size)? new_ranks[i+1] : MPI_PROC_NULL;
    printf("[%d] %d p %d n %d", world_rank, i, previous, next);

    if ( world_rank == 0 ) {
        /* 3. Etiquetado inicial */
        matrixResult= (int *)malloc( (rows)*(columns) * sizeof(int) );
        matrixResultCopy= (int *)malloc( (rows)*(columns) * sizeof(int) );
        if ( (matrixResult == NULL)  || (matrixResultCopy == NULL)  ) {
            perror ("Error reservando memoria");
            return -1;
        }
        /*for(i=0;i< rows; i++){
          for(j=0;j< columns; j++){
          matrixResult[i*(columns)+j]=-1;
        // Si es 0 se trata del fondo y no lo computamos
        if(matrixData[i*(columns)+j]!=0){
        matrixResult[i*(columns)+j]=i*(columns)+j;
        }
        }
        }*/
        // Pack dimensions of matrix
        dimensions[0] = rows;
        dimensions[1] = columns;


#ifdef DEBUG
        printf("Data\n");
        for(i=0; i<rows; ++i){
            for(j=0; j<columns; ++j){
                printf("%d ", matrixData[i*columns+j]);
            }
            printf("\n");
        }
        fflush(stdout);
        printf("Matrix\n");
        /*for(i=0; i<rows; ++i){
          for(j=0; j<columns; ++j){
          printf("%d ", matrixResult[i*columns+j]);
          }
          printf("\n");
          }
          fflush(stdout);*/
#endif
    }


    // Broadcast matrix dimensions
    MPI_Bcast(dimensions, 2, MPI_INT, 0, MPI_COMM_WORLD);

    // Update matrix rows for each process
    if( world_rank ){
        rows = dimensions[0];
        columns = dimensions[1];
    }

    if( world_rank == 0){
        // FIXME: Incorrect way to calculate displacements
        // Calculate vector of rows for each proc with borders

        for(i=0; i<world_size; i++){
            j = new_ranks[i];
            vectorRows[j] = rows/world_size+(i < rows%world_size? 1: 0);
        }
        /*
           vectorRows[i] = rows/world_size+(rows*world_size%world_size);
           vectorRows[world_size-1] += rows - ( rows/world_size+(rows*world_size%world_size))*world_size;
         */
        // Calculate number of cells for each process
        for(i=0; i<world_size; i++){
            vectorSizes[i] = vectorRows[i]*columns;
        }

        // Calculate displacements for each process
        vectorDis[0] = 0;
        for(i=1; i<world_size; i++){
            j = new_ranks[i];
            vectorDis[j] = vectorDis[new_ranks[i-1]]+vectorSizes[new_ranks[i-1]];
        }

#ifdef DEBUG
        for(i=0; i<world_size; i++)
            printf("#%d[%d] r %d d %d s %d",i, new_ranks[i], vectorRows[new_ranks[i]], vectorDis[new_ranks[i]], vectorSizes[new_ranks[i]]);
        printf("\n");fflush(stdout);
#endif

        if( world_size > 1){
            // Update vector of rows for each proc with borders
            vectorRows[0] += 1;
            vectorRows[new_ranks[world_size-1]] += 1;
            for(i=1; i<world_size-1; i++){
                j = new_ranks[i];
                vectorRows[j] += 2;
            }

        }

        // Update number of cells for each process
        for(i=0; i<world_size; i++){
            j = new_ranks[i];
            vectorSizes[i] = vectorRows[i]*columns;
        }

    }

    // TODO: Avoid this communication. Proposal:
    // 1. Each proc calculate its number of rows(`sub_rows`) by knowing `rows` and `world_rank`
    // 2. Allgather of `sub_rows` into `vectorRows`, which its in all procs but only used in rank=0.
    // 3. Each proc calculate `displacement` knowing `sub_rows` of procs from rank=0 to rank=previous.
    //    3.1. rank=0 stays calculating vectorSizes and vectorDis
    // 4. Avoid Scatter of `displacement` as each proc already knows it.
    // 5. Each proc update number of `sub_rows` with borders.
    // 6. Each proc calculate ncells, which has to be equal as `vectorSizes[rank]` at rank=0.

    // Scatter number of rows for each proc
    MPI_Scatter(vectorRows, 1, MPI_INT, &sub_rows, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate number of cells for each submatrix
    ncells = sub_rows*columns;

    if(world_size > 1){
        // Assign MPI_PROC_NULL wildcards for communications at first and last proc
        if( world_rank == 0 ){
            previous = MPI_PROC_NULL;
        } else if( world_rank == world_size-1) {
            next = MPI_PROC_NULL; 
        }
    }

    // Scatter each displacement
    MPI_Scatter(vectorDis, 1, MPI_INT, &displacement, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(world_rank==0){
        // Move displacement one row upwards for each proc
        for(i=1; i<world_size; i++)
            vectorDis[i] -= columns;
#ifdef DEBUG
        for(i=0; i<world_size; i++)
            printf("#%d[%d] r %d d %d s %d",i, new_ranks[i], vectorRows[new_ranks[i]], vectorDis[new_ranks[i]], vectorSizes[new_ranks[i]]);
        printf("\n");fflush(stdout);
#endif
    }

    sub_matrixData = (int*) malloc(sizeof(int)*sub_rows*columns);
    sub_matrixResult = (int*) malloc(sizeof(int)*sub_rows*columns);
    sub_matrixResultCopy = (int*) malloc(sizeof(int)*sub_rows*columns);

    // Broadcast of data matrix
    // MPI_Bcast(matrixData, rows*columns, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatterv(matrixData, vectorSizes, vectorDis, MPI_INT, sub_matrixData, ncells, MPI_INT, 0, MPI_COMM_WORLD);
    
    // FIXME: Unused
    //MPI_Scatterv(matrixResult, vectorSizes, vectorDis, MPI_INT, sub_matrixResult, ncells, MPI_INT, 0, MPI_COMM_WORLD);

    // Fill submatrix for ResultCopy with -1
    for(i=0; i<sub_rows; ++i){
        for(j=0; j<columns; ++j){
            sub_matrixResultCopy[i*columns+j]=-1;
            sub_matrixResult[i*(columns)+j]=-1;
            // Si es 0 se trata del fondo y no lo computamos
            if(sub_matrixData[i*(columns)+j]!=0){
                sub_matrixResult[i*(columns)+j]=(world_rank==0?i:i-1)*columns+j+displacement;
            }
        }
    }

#ifdef DEBUG
    if(world_rank==0){
        printf("[%d] Init submatrix of Data\n", world_rank);
        for(i=0; i<sub_rows; ++i){
            for(j=0; j<columns; ++j){
                printf("%d ", sub_matrixData[i*columns+j]);
            }
            printf("\n");
        }
        fflush(stdout);
        printf("[%d] Init submatrix of Result\n", world_rank);
        for(i=0; i<sub_rows; ++i){
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
    int flagCambio=1, flagCambioGlobal=1; 

    /* 4.2 Busqueda de los bloques similiares */
    for(t=0; flagCambioGlobal !=0; t++){
        printf("Init %d #%d\n", world_rank, t);

        if(world_size > 1){
            printf("Sending %d #%d %d n %d\n", world_rank, t, previous, next);
            // Segment fault at sending on both programs
            MPI_Isend(sub_matrixResult+columns, columns, MPI_INT, previous, flagCambio, MPI_COMM_WORLD, &reqs[0]);  // Send row #1 to previous proc
            MPI_Isend(sub_matrixResult+(sub_rows-2)*columns, columns, MPI_INT, next, flagCambio, MPI_COMM_WORLD, &reqs[1]);  // Send penultimate row to next proc
            printf("Receiving %d #%d %d n %d\n", world_rank, t, previous, next);
            MPI_Irecv(sub_matrixResultCopy+(sub_rows-1)*columns, columns, MPI_INT, next, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[0]);  // Wait to last row from next proc
            MPI_Irecv(sub_matrixResultCopy, columns, MPI_INT, previous, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[1]);  // Wait to first row from previous pro

#ifdef DEBUG
            if(world_rank==0){
                printf("[%d]%d Submatrix\n",world_rank, t);
                for(i=0; i<sub_rows; ++i){
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
        for(i=1;i<sub_rows-1;i++){
            for(j=1;j<columns-1;j++){
                if(sub_matrixResult[i*(columns)+j]!=-1){
                    sub_matrixResultCopy[i*(columns)+j]=sub_matrixResult[i*(columns)+j];
                }
            }
        }
        printf("Waiting %d #%d to %d n %d\n", world_rank, t, previous, next);
        if(world_size > 1){
            MPI_Waitall(2, reqs, stats);
        }
#ifdef DEBUG
        if(world_rank==0){
            printf("[%d]%d Submatrix copy\n",world_rank, t);
            for(i=0; i<sub_rows; ++i){
                for(j=0; j<columns; ++j){
                    printf("%d ", sub_matrixResultCopy[i*columns+j]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
#endif

        printf("Let's go %d #%d\n", world_rank, t);
        if(flagCambio || 
                (stats[0].MPI_TAG && next!=MPI_PROC_NULL) || 
                (stats[1].MPI_TAG && previous!=MPI_PROC_NULL)
                ){
            flagCambio=0; 
            /* 4.2.2 Computo y detecto si ha habido cambios */
            for(i=1;i<sub_rows-1;i++){
                for(k=1;k<columns-1;k++){
                    int result,sol;
                    j=i*columns+k;
                    // int dataj = (world_rank==0?i:i-1)*columns+k+displacement;
                    result = sub_matrixResultCopy[j];
                    sol=0;
                    if(result!=-1){
                        //Si es de mi mismo grupo, entonces actualizo
                        if(sub_matrixData[j-columns] == sub_matrixData[j])
                        {
                            result = min (result, sub_matrixResultCopy[j-columns]);
                        }
                        if(sub_matrixData[j+columns] == sub_matrixData[j])
                        {
                            result = min (result, sub_matrixResultCopy[j+columns]);
                        }
                        if(sub_matrixData[j-1] == sub_matrixData[j])
                        {
                            result = min (result, sub_matrixResultCopy[j-1]);
                        }
                        if(sub_matrixData[j+1] == sub_matrixData[j])
                        {
                            result = min (result, sub_matrixResultCopy[j+1]);
                        }
                        // Si el indice no ha cambiado retorna 0
                        if(sub_matrixResult[j] == result){ sol=0; }
                        // Si el indice cambia, actualizo matrixResult de resultados con el indice adecuado y retorno 1
                        else { sub_matrixResult[j]=result; sol=1;}
                        flagCambio|=sol;
                    }
                }
            }
        }

#ifdef DEBUG
        if(world_rank==0){
            printf("[%d]%d Submatrix result\n",world_rank, t);
            for(i=0; i<sub_rows; ++i){
                for(j=0; j<columns; ++j){
                    printf("%d ", sub_matrixResult[i*columns+j]);
                }
                printf("\n");
            }
            fflush(stdout);
        }
#endif

        MPI_Allreduce(&flagCambio, &flagCambioGlobal, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
        printf("Finish %d #%d\n", world_rank, t);
    }

    // MPI_Gatherv(sub_matrixResult+(world_rank?columns:0), ncells, MPI_INT, matrixResult, vectorSizes, vectorDis, MPI_INT, 0, MPI_COMM_WORLD);   

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
    for(i=1;i<sub_rows-1;i++){
        for(j=1;j<columns-1;j++){
            if(sub_matrixResult[i*columns+j] == (world_rank==0?i:i-1)*columns+j+displacement) nBlocks++; 
        }
    }        

    MPI_Reduce(&nBlocks, &numBlocks, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

#ifdef DEBUG
    // Substract borders from number of cells
    if(world_rank==0||world_rank==world_size-1)
        ncells -= columns;
    else
        ncells -= columns*2;
    // Gather number of cells and displacement of each proc without borders
    MPI_Gather(&ncells, 1, MPI_INT, vectorSizes, 1, MPI_INT, 0, MPI_COMM_WORLD); 
    MPI_Gather(&displacement, 1, MPI_INT, vectorDis, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    // Gather result matrix
    MPI_Gatherv(sub_matrixResult+(world_rank==0? 0: columns), ncells, MPI_INT, matrixResult, vectorSizes, vectorDis, MPI_INT, 0, MPI_COMM_WORLD);   
    if(world_rank==0){
        printf("Matrix result\n");
        for(i=0;i<rows;i++){
            for(j=0;j<columns;j++){
                printf ("%d\t", matrixResult[i*columns+j]);
            }
            printf("\n");
        }
    }
#endif

#ifdef DEBUG
    printf("[%d]subs dataX%08x resultX%08x resultCX%08x\n", world_rank, sub_matrixData, sub_matrixResult, sub_matrixResultCopy);fflush(stdout);
#endif
    free(sub_matrixData);
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
