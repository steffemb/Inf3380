#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>


void print_mat (double*** matrix, int* num_rows, int* num_cols);
void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols);
void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* num_rows, int* num_cols);
void write_matrix_binaryformat (char* filename, double** matrix,
                                int num_rows, int num_cols);
void allocate_matrix(double*** matrix, int* num_rows, int* num_cols);
int free2dfloat(double ***array);


int main(void)
{

    clock_t start = clock(), diff;
    double **matrix_a;
    double **matrix_b;
    double **matrix_c;
    int num_rows = 100; // for a, switched for b.
    int num_cols = 50;
    int i,j,k;
    double dummy = 0;
    //int counter = 0;

    int th_id, num_threads;

    read_matrix_binaryformat ("small_matrix_a.bin", &matrix_a, &num_rows, &num_cols); //get matrix a
    read_matrix_binaryformat ("small_matrix_b.bin", &matrix_b, &num_cols, &num_rows); //get matrix b
    allocate_matrix(&matrix_c, &num_rows, &num_rows);

    //#pragma omp parallel num_threads (4) private(matrix_a) private(matrix_b) shared(matrix_c)

//#pragma omp parallel default(shared) private(th_id,i,j,k,dummy)
//    {
//        th_id = omp_get_thread_num();
//#pragma omp barrier
//        if ( th_id == 0 ) {
//            num_threads = omp_get_num_threads();
//            printf("There are %d threads\n",num_threads);

 //       }
#pragma omp parallel for num_threads(2) default(shared) private(i,j,k)

        for (i=0; i<(num_rows); i++){

            for (j=0; j<(num_rows); j++){
                for (k=0; k<(num_cols); k++){
                    (matrix_c)[i][j] += (matrix_a)[i][k]*(matrix_b)[k][j];

                }//k

                //printf(" dummy = %f \n", dummy);

                //counter += 1;
            }//j
            //printf("%d \n",counter);
        }//i


    //matrix_multiplication(&matrix_a, &matrix_b, &matrix_c, &num_rows, &num_cols); // compute c

    //print_mat (&matrix_a, &num_rows, &num_cols);
    //print_mat (&matrix_b, &num_rows, &num_cols);
    //print_mat (&matrix_c, &num_rows, &num_rows);
     write_matrix_binaryformat ("my_matrix_c.bin", matrix_c, num_rows, num_rows);
    diff = clock() - start;
    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken %d seconds %d milliseconds \n", msec/1000, msec%1000);

    return 0;
}

void matrix_multiplication(double*** matrix_a, double*** matrix_b,
                           double*** matrix_c, int* num_rows, int* num_cols){
    int i,j,k;
    double dummy = 0;
    int counter = 0;

    /* storage allocation of matrix c */
    *matrix_c = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix_c)[0] = (double*)malloc((*num_rows)*(*num_rows)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix_c)[i] = (*matrix_c)[i-1]+(*num_rows);


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
        for (j=1; j<(*num_cols); j++){
            printf("%f ", (*matrix)[i][j]); // i = 100, j = 50
            //printf("\n");
        }
    }
}


void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols)
{
    int i;
    FILE* fp = fopen (filename,"rb");
    fread (num_rows, sizeof(int), 1, fp);
    fread (num_cols, sizeof(int), 1, fp);

    /* storage allocation of the matrix */
    *matrix = (double**)malloc((*num_rows)*sizeof(double*));
    (*matrix)[0] = (double*)malloc((*num_rows)*(*num_cols)*sizeof(double));
    for (i=1; i<(*num_rows); i++)
        (*matrix)[i] = (*matrix)[i-1]+(*num_cols);

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

int free2dfloat(double ***array) {
    /* free the memory - the first element of the array is at the start */
    free(&((*array)[0][0]));

    /* free the pointers into the memory */
    free(*array);

    return 0;
}

