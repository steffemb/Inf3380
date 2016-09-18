#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>


void print_mat (double*** matrix, int* num_rows, int* num_cols);
void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols);
void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* a_num_rows, int* b_num_rows, int *b_num_cols);
void write_matrix_binaryformat (char* filename, double** matrix,
                                int num_rows, int num_cols);
void allocate_matrix(double*** matrix, int* num_rows, int* num_cols);
void allocate_vector(double** vector, int* length);

int free2dfloat(double ***array);

void sammenlikn(double*** array_a, double*** array_b, int* num_rows, int* num_cols, int* value);


int main (int argc, char *argv[])
{
    clock_t start = clock(), diff;
    int error = 0;
    int print_state = 0;
    char *input_matrix_a, *input_matrix_b;
    double **matrix_a;
    double **matrix_b;
    double **matrix_c;
    double **my_matrix_b;
    double **my_matrix_c;
    double ** test_matrix_c;
    int num_rows = 100; // for a, switched for b.
    int num_cols = 50;
    int my_num_rows_b, my_num_cols_b, my_num_rows_c, my_num_cols_c;
    int mat_size = num_rows*num_cols;
    int my_rank, num_procs, k, i, j, input_state;
    int value = 0;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //input handling

    if (argc == 2){
        if (!strcmp(argv[1],"-help")){
            input_state = 0;
            error = 1;
        }else if (!strcmp(argv[1],"-print_matrix")){
            input_state = 0;
            print_state = 1;
            input_state = 1;
            input_matrix_a = "small_matrix_a.bin";
            input_matrix_b = "small_matrix_b.bin";

        }else if (!strcmp(argv[1],"-bench")){
            input_state = 0;
            print_state = 2;
            input_state = 1;
            input_matrix_a = "large_matrix_a.bin";
            input_matrix_b = "large_matrix_b.bin";
            num_rows = 1000;
            num_cols = 500;
        }else{
            input_state = 2;
            error = 1;
        }

    }else{
        if (argc == 1){//defaults
            input_state = 1;
            input_matrix_a = "small_matrix_a.bin";
            input_matrix_b = "small_matrix_b.bin";
        }
    }
    if (my_rank==0){ //printouts of input params and state
        if(input_state == 0){
            printf("-----------------------help---------------------\n");
            printf("\n");
            printf("program takes 0,1 or two arguments. 0-->default, 1: -print_matrix or -bench");
            printf("the last option is to specify input and output filenames(not implemented yet)");
            printf("\n");
        }
        if(input_state == 1){
            printf("\n");
            printf("Using default values, small_matrix_a.bin and small_matrix_a.bin");
            printf("\n");
        }
        if(input_state == 2){
            printf("\n");
            printf("Error: try -help --EXITING--\n");
            printf("\n");
        }
    }

    //take care of invalid input and abort
    if (error != 0) {
        if (my_rank == 0) {
            fprintf(stderr, "Error: Program terminated with error code %d\n", error);
        }
        MPI_Finalize ();
        exit(error);
    }

    //end of input handling


    //printf("my rank = %d, of %d \n", my_rank, num_procs); //report

    my_num_rows_b = num_cols; //block partition only on cols
    my_num_cols_b = num_rows/(num_procs);
    my_num_rows_c = num_rows; //block partition only on cols
    my_num_cols_c = num_rows/(num_procs);

    allocate_matrix(&matrix_a, &num_rows, &num_cols); //every process will need matrix a, therefore outside the "if" (can be done more memory effecticve?)
    allocate_matrix(&my_matrix_c, &my_num_rows_c, &my_num_cols_c);
    allocate_matrix(&my_matrix_b, &my_num_rows_b, &my_num_cols_b); //give local b same dimensions as local a. This simplifies all-to-all


    //let one process read files
    if(my_rank == 0){
        allocate_matrix(&matrix_b, &num_cols, &num_rows);
        allocate_matrix(&matrix_c, &num_rows, &num_rows);

        read_matrix_binaryformat (input_matrix_a, &matrix_a, &num_rows, &num_cols); //get matrix a
        read_matrix_binaryformat (input_matrix_b, &matrix_b, &num_cols, &num_rows); //get matrix b (split these on two procs?)

    }



    MPI_Bcast(&(matrix_a[0][0]),mat_size,MPI_DOUBLE,0,MPI_COMM_WORLD); //every proc has a


    //-------------------------------------------------
    //for use with send - recv, will try scatter later
    if(my_rank != 0){
        for(i = 0; i< my_num_rows_b; i++){
            MPI_Recv(& my_matrix_b[i][0], my_num_cols_b , MPI_DOUBLE,0, 100,MPI_COMM_WORLD, &status);

        }
    }
    if (my_rank == 0){

        //keep my part of b

        for(k = 0; k< my_num_rows_b; k++){
            for(j = 0; j< my_num_cols_b; j++){
                my_matrix_b[k][j] = matrix_b[k][j];
            }
            for(i = 1; i< num_procs; i++){
                MPI_Send(& matrix_b[k][i*my_num_cols_b], my_num_cols_b , MPI_DOUBLE, i, 100, MPI_COMM_WORLD);
            }
        }
    }
    //-------------------------------

    //MPI_Scatter((matrix_b), my_num_rows_b, MPI_DOUBLE, (my_matrix_b), my_num_rows_b, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //for(k=0; k<my_num_rows_b; k++){
    //    MPI_Scatter(&(matrix_b[k][0]), my_num_cols_b, MPI_DOUBLE, &(my_matrix_b[k][0]), my_num_cols_b, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //}


    matrix_multiplication(&matrix_a, &my_matrix_b, &my_matrix_c, &num_rows, &my_num_rows_b, &my_num_cols_b);
    //print_mat (&my_matrix_c, &num_rows, &my_num_cols_b);

    if(my_rank != 0){

        for(i = 0; i< num_rows; i++){
            MPI_Send(& my_matrix_c[i][0], my_num_cols_b , MPI_DOUBLE, 0, 100, MPI_COMM_WORLD);
        }
    }

    if (my_rank == 0){
        for(k = 0; k< num_rows; k++){
            for(j = 0; j< my_num_cols_b; j++){
                matrix_c[k][j] = my_matrix_c[k][j];
            }
            for(i = 1; i< num_procs; i++){
                MPI_Recv(& matrix_c[k][i*my_num_cols_b], my_num_cols_b , MPI_DOUBLE, i, 100, MPI_COMM_WORLD, &status);
            }
        }
    }


    if(my_rank == 0){
        if(print_state == 1){
            //print and check if computation is correct
            print_mat (&matrix_c, &num_rows, &num_rows);
            allocate_matrix(&test_matrix_c, &num_rows, &num_rows);
            read_matrix_binaryformat ("small_matrix_c.bin", &test_matrix_c, &num_rows, &num_cols);
            sammenlikn(&matrix_c, &test_matrix_c, &num_rows, &num_cols,&value);
            free2dfloat(&test_matrix_c);
        }
        //print_mat (&my_matrix_b, &my_num_rows_b, &my_num_cols_b);
        write_matrix_binaryformat ("my_matrix_c.bin", matrix_c, num_rows, num_rows);
    }

    MPI_Barrier (MPI_COMM_WORLD); //no free 4 u before it used!

    free2dfloat(&matrix_a);
    free2dfloat(&my_matrix_c);
    free2dfloat(&my_matrix_b);
    if(my_rank == 0){
        free2dfloat(&matrix_b);
        free2dfloat(&matrix_c);
    }

    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    if(print_state == 2){
        printf("Time taken %d seconds %d milliseconds, I am proc %d \n", msec/1000, msec%1000, my_rank);
    }
    MPI_Finalize ();
    return error;
}

void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* a_num_rows, int* b_num_rows, int* b_num_cols){
    int i,j,k;
    double dummy = 0;
    int counter = 0;

    for (i=0; i<(*a_num_rows); i++){//100
        for (j=0; j<(*b_num_cols); j++){//25
            for (k=0; k<(*b_num_rows); k++){//50
                dummy += (*matrix_a)[i][k]*(*matrix_b)[k][j];

            }//k
            (*matrix_c)[i][j] = dummy;
            //printf(" dummy = %f \n", dummy);
            dummy = 0;
            counter += 1;
        }//j
        //printf("%d \n",counter);
    }//i
}

void write_matrix_binaryformat (char* filename, double** matrix,
                                int num_rows, int num_cols)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&num_rows, sizeof(int), 1, fp);
    fwrite (&num_cols, sizeof(int), 1, fp);
    fwrite (matrix[0], sizeof(double), (num_rows)*(num_cols), fp);
    fclose (fp);
}

void print_mat (double*** matrix, int* num_rows, int* num_cols){
    printf("%d \n", *num_rows);
    printf("%d \n", *num_cols);


    int i, j;
    for (i=1; i<(*num_rows); i++){
        printf("\n");
        printf("\n");
        for (j=1; j<(*num_cols); j++){
            printf("%f ", (*matrix)[i][j]); // i = 100, j = 50
            //printf("\n");
        }
    }
}


void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols)
{
    FILE* fp = fopen (filename,"rb");
    fread (num_rows, sizeof(int), 1, fp);
    fread (num_cols, sizeof(int), 1, fp);



    /* read in the entire matrix */
    fread ((*matrix)[0], sizeof(double), (*num_rows)*(*num_cols), fp);
    fclose (fp);
}

void allocate_vector(double** vector, int* length){
    int i;
    /* storage allocation of the matrix */
    *vector = (double*)malloc((*length)*sizeof(double));

}



void allocate_matrix(double*** matrix, int* num_rows, int* num_cols){
    int i;
    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
}

int free2dfloat(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

void sammenlikn(double*** array_a, double*** array_b, int* num_rows, int* num_cols,int* value){
    int i,j;

    int counter = 0;
    for (i=0; i<(*num_rows); i++){

        for (j=0; j<(*num_cols); j++){

            if ( (*array_a)[i][j] =! (*array_b)[i][j] ){
                counter += 1;
                printf("error");
            }
        }
    }
    *value = counter;
    printf("\n");
    printf("\n");
    printf("number of errors = %d \n", counter);
    printf("\n");
}
