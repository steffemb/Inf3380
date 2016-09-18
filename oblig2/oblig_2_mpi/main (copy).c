#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>


void print_mat (double*** matrix, int* num_rows, int* num_cols);
void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols);
void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* num_rows, int* num_cols);
void write_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols);
void allocate_matrix(double*** matrix, int* num_rows, int* num_cols);


int main (int argc, char *argv[])
{
    double **matrix_a;
    double **matrix_b;
    double **matrix_c;
    int num_rows = 100; // for a, switched for b.
    int num_cols = 50;
    int mat_size = num_rows*num_cols;
    int my_rank, num_procs;

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    printf("my rank = %d, of %d \n", my_rank, num_procs);

    allocate_matrix(&matrix_a, &num_rows, &num_cols);
    allocate_matrix(&matrix_b, &num_cols, &num_rows);
    allocate_matrix(&matrix_c, &num_rows, &num_rows);

    //let one process read files
    if(my_rank == 0){
        read_matrix_binaryformat ("small_matrix_a.bin", &matrix_a, &num_rows, &num_cols); //get matrix a
        read_matrix_binaryformat ("small_matrix_b.bin", &matrix_b, &num_cols, &num_rows); //get matrix b

    }
    //for now we will use 100 processes, which get one row of a, and the whole of b
    MPI_Bcast(&matrix_b[0][0],mat_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&matrix_a[0][0],mat_size, MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Barrier (MPI_COMM_WORLD);


    //matrix_multiplication(&matrix_a, &matrix_b, &matrix_c, &num_rows, &num_cols); // compute c

    //print_mat (&matrix_a, &num_rows, &num_cols);
    //print_mat (&matrix_b, &num_rows, &num_cols);
    if(my_rank == 1){

        matrix_multiplication(&matrix_a, &matrix_b, &matrix_c, &num_rows, &num_cols); // compute c
        print_mat (&matrix_c, &num_rows, &num_rows);
        write_matrix_binaryformat ("small_matrix_c.bin", &matrix_c, &num_rows, &num_rows);
    }

    MPI_Barrier (MPI_COMM_WORLD);
    MPI_Finalize ();
    return 0;
}

void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* num_rows, int* num_cols){
    int i,j,k;
    double dummy = 0;
    int counter = 0;



    for (i=0; i<(*num_rows); i++){

        for (j=0; j<(*num_rows); j++){
            for (k=0; k<(*num_cols); k++){
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

void write_matrix_binaryformat (char* filename, double*** matrix,
                                int* num_rows, int* num_cols)
{
    FILE *fp = fopen (filename,"wb");
    fwrite (&num_rows, sizeof(int), 1, fp);
    fwrite (&num_cols, sizeof(int), 1, fp);
    fwrite (&matrix[0], sizeof(double), (*num_rows)*(*num_cols), fp);
    fclose (fp);
}

void print_mat (double*** matrix, int* num_rows, int* num_cols){
    printf("%d \n", *num_rows);
    printf("%d \n", *num_cols);


    int i, j;
    for (i=1; i<(*num_rows); i++){
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

void allocate_matrix(double*** matrix, int* num_rows, int* num_cols){
    int i;
    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(*num_cols);
}
