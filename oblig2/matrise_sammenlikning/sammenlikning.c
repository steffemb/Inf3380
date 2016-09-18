#include <stdio.h>
#include <stdlib.h>



void print_mat (double*** matrix, int* num_rows, int* num_cols);
void read_matrix_binaryformat (char* filename, double*** matrix,
                               int* num_rows, int* num_cols);
void allocate_matrix(double*** matrix, int* num_rows, int* num_cols);
int free2dfloat(double ***array);
void sammenlikn(double*** array_a, double*** array_b, int* num_rows, int* num_cols, int* value);


int main (int argc, char *argv[]){
    double **matrix_a;
    double **matrix_b;
    int num_rows = 100;
    int num_cols = 100;
    int value = 0;
    char *input_matrix_a, *input_matrix_b;

    input_matrix_a = "matrix_c.bin";
    input_matrix_b = "my_matrix_c.bin";


    allocate_matrix(&matrix_a, &num_rows, &num_cols);
    allocate_matrix(&matrix_b, &num_rows, &num_cols);

    read_matrix_binaryformat (input_matrix_a, &matrix_a, &num_rows, &num_cols);
    read_matrix_binaryformat (input_matrix_b, &matrix_b, &num_rows, &num_cols);

    print_mat (&matrix_b, &num_rows, &num_rows);

    sammenlikn(&matrix_a, &matrix_b, &num_rows, &num_cols, &value);




    return 0;
}

void sammenlikn(double*** array_a, double*** array_b, int* num_rows, int* num_cols,int* value){
    int i,j;
    double eps = 0.001;
    int counter = 0;
    for (i=0; i<(*num_rows); i++){

        for (j=0; j<(*num_cols); j++){


            if ( (*array_a)[i][j] < (*array_b)[i][j] + eps && (*array_a)[i][j] + eps > (*array_b)[i][j] ){
                counter += 1;
                //printf("error\n");
            }
        }
    }
    printf("\n");
    printf("number of errors = %d\n", counter);
    *value = counter;
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
