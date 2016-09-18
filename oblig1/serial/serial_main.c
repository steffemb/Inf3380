#include "../simple-jpeg/import_export_jpeg.h"
#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>



typedef struct{
    float ** image_data;
    int m, n;
}image; // define how variable type "image" should be

void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char* image_chars, image *u);
void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);


int main (int argc, char *argv[]) {
    char *input_jpeg_filename, *output_jpeg_filename;
    int image_height, image_width, num_components;
    int iters;// = atoi(argv[2]);
    float kappa;// = atof(argv[1]);
    unsigned char* image_chars;
    image u, u_bar; //typedef struct


    //handling input params, take care of any posible input

    if (argc == 2){
        if (!strcmp(argv[1],"help")){
            printf("-----------------------help---------------------\n");
            printf("\n");
            printf("Program for applying ISO-Diffusion to JPEG\n");
            printf("usage: execute with two or four arguments.\n");
            printf("for two args, they must be iterations + kappa.\n");
            printf("for four args, add file_in and file_out.\n");
            printf("ex: computer@username: serial_main 20 1.0 file_in.jpeg file_out.jpeg\n");
            printf("\n");
            exit(EXIT_SUCCESS);
        }else{
            printf("Error: give two, or four arguments, order: kappa, iters, innput_file, output_file. --EXITING--\n");
            exit(EXIT_SUCCESS);
        }
    }else if (argc == 1){//defaults
        printf("Using default values, kappa = 0.1, iters = 10\n");
        kappa = 0.1;
        iters = 10;
        input_jpeg_filename = "mona_lisa_noisy.jpeg";
        output_jpeg_filename = "mona_lisa_smothed.jpeg";
    }else if(argc ==3){
        if (atof(argv[2])>0 && atof(argv[2])<1) {
            printf("Using default filenames\n");
            int iters = atoi(argv[1]);
            float kappa = atof(argv[2]);
            input_jpeg_filename = "mona_lisa_noisy.jpeg";
            output_jpeg_filename = "mona_lisa_smothed.jpeg";
        }else{
            printf("kappa must be in the range 0-1, --EXITING--\n");
            exit(EXIT_SUCCESS);
        }
    }else if (argc == 4){
        printf("Error: give two, or four arguments, order: kappa, iters, innput_file, output_file. --EXITING--\n");
        exit(EXIT_SUCCESS);
    }else if (argc == 5 ){
        if (atof(argv[2])>0 && atof(argv[2])<1) {
            input_jpeg_filename = argv[3];
            output_jpeg_filename = argv[4];
            int iters = atoi(argv[1]);
            float kappa = atof(argv[2]);
            printf("set: kappa= %f, iterations =  %d, input= %s, output = %s\n", atof(argv[1]), atoi(argv[2]), argv[3], argv[4]);
        }else{
            printf("kappa must be in the range 0-1, --EXITING--\n");
            exit(EXIT_SUCCESS);
        }
    }else if (argc > 5 ){
        printf("too many arguments. max = 4, --EXITING--\n");
        exit(EXIT_SUCCESS);
    }

    //end of input handling

    // get 1d data
    import_JPEG_file(input_jpeg_filename, &image_chars, &image_height, &image_width, &num_components);

    //allocate memory and fill with data
    allocate_image(&u, image_height, image_width);
    allocate_image(&u_bar, image_height, image_width);
    //1d ---> 2d
    convert_jpeg_to_image (image_chars, &u);
    //printf("W:%d, H:%d, C:%d name:%s ",image_width, image_height, num_components, input_jpeg_filename );

    //smoooothen
    iso_diffusion_denoising (&u, &u_bar, kappa, iters);

    //2d ---> 1d
    convert_image_to_jpeg (&u_bar, image_chars);
    export_JPEG_file(output_jpeg_filename, image_chars, image_height, image_width, num_components, 75);

    //free memory
    deallocate_image(&u);
    deallocate_image (&u_bar);

    return 0;
}

//functions

void allocate_image(image *u, int m, int n){
    /* allocates memory for a struct image */
    int i;
    u->m = m;
    u->n = n;
    u->image_data = (float**)malloc(m*sizeof(float*));
    for (i=0; i<m; i++){
        u->image_data[i] = (float*)malloc(n*sizeof(float*));
    }
}


void deallocate_image(image *u){
    /* deallocates memory of a struct image */
    int i;
    for (i=0; i<u->m; i++){
        free(u->image_data[i]);
    }
    free(u->image_data);
}


void convert_jpeg_to_image(const unsigned char* image_chars, image *u){
    /* converts an array of *constant char data to struct image data */
    int i, j, k;
    k = 0;
    for (i=0; i<u->m; i++){
        for (j=0; j<u->n; j++){
            // is i_max = 4289, and j_max = 2835???
            //j = floor((i/(u->n-1)));
            u->image_data[i][j] = ((float)image_chars[k])/255; // woooorks??? [j][i]??
            k += 1; //check k count if jibberjabber
        }
    }
}




void convert_image_to_jpeg(const image *u, unsigned char* image_chars){
    int i, j, k;
    k = 0;
    for (i=0; i<u->m; i++){
        for (j=0; j<u->n; j++){
            image_chars[k] =  (unsigned char)((int)(u->image_data[i][j]*255));
            k += 1;
        }
    }
}



void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters){
    /* applies iso diffusion to a struct image (2d array, or "matrix") */
    /* this function leaves the boundary values untouched, but we could */
    /* have implemented a 3-point iso diffusion explicitly on the boundarys */
    /* but this woulg give allot of overhead */

    int i, j, k;
    for (i=0; i<u->m; i++){
        for (j=0; j<u->n; j++){ //check max vals, is i_max = 4289, and j_max = 2835???
            u_bar->image_data[i][j] = u->image_data[i][j];
        }
    }
    for (k=0; k<iters; k++){
        for (i=1; i<u->m-1; i++){
            for (j=1; j<u->n-1; j++){
                u_bar->image_data[i][j] = (u->image_data[i][j] + kappa*(u->image_data[i-1][j] + u->image_data[i][j-1] - (4*u->image_data[i][j]) + u->image_data[i][j+1] + u->image_data[i+1][j]));

            }
        }

        for (i=0; i<u->m; i++){
            for (j=0; j<u->n; j++){// --copy for new iter
                u->image_data[i][j] = u_bar->image_data[i][j];


            }
        }
    }//k loop
}
