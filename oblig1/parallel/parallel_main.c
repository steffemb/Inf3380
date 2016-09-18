
#include "../simple-jpeg/import_export_jpeg.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>



typedef struct{
    float ** image_data;
    int m, n;
}image; // define how variable type "image" should be

void allocate_image(image *u, int m, int n);
void deallocate_image(image *u);
void convert_jpeg_to_image(const unsigned char *image_chars, image *u);
void iso_diffusion_denoising (image *u, image *u_bar, float kappa, int iters);
void convert_image_to_jpeg(const image *u, unsigned char* image_chars);


int main (int argc, char *argv[]) {

    int error = 0;
    int input_state, dummy;
    char *input_jpeg_filename, *output_jpeg_filename;
    int image_height, image_width, num_components, iters;
    int my_m, my_n, my_rank, num_procs, i, k;
    float kappa;
    unsigned char *image_chars;
    image u, u_bar, whole_image; //typedef struct

    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    //printf("my rank is %d \n", my_rank);

    //handling input params, take care of any posible input
    if (argc == 2){
        if (!strcmp(argv[1],"help")){
            input_state = 0;
            error = 1;
        }else{
            input_state = 2;
            error = 1;
        }

    }else{
        if (argc == 1){//defaults
            input_state = 1;
            kappa = 0.1;
            iters = 10;
            input_jpeg_filename = "mona_lisa_noisy.jpeg";
            output_jpeg_filename = "mona_lisa_smothed.jpeg";
        }else if(argc ==3){
            if (atof(argv[2])>0 && atof(argv[2])<1) {
                input_state = 3;
                iters = atoi(argv[1]);
                kappa = atof(argv[2]);
                input_jpeg_filename = "mona_lisa_noisy.jpeg";
                output_jpeg_filename = "mona_lisa_smothed.jpeg";
            }else{
                input_state = 4;
                error = 1;
            }
        }else if (argc == 4){
            input_state = 5;
            error = 1;
        }else if (argc == 5 ){
            if (atof(argv[2])>0 && atof(argv[2])<1) {
                input_state = 6;
                input_jpeg_filename = argv[3];
                output_jpeg_filename = argv[4];
                iters = atoi(argv[1]);
                kappa = atof(argv[2]);

            }else{
                input_state = 7;
                error = 1;
            }
        }else if (argc > 5 ){
            input_state = 8;
            error = 1;
        }
    }
    if (my_rank==0){ //printouts of input params and state
        if(input_state == 0){
            printf("-----------------------help---------------------\n");
            printf("\n");
            printf("Program for applying ISO-Diffusion to JPEG\n");
            printf("usage: execute with two or four arguments.\n");
            printf("for two args, they must be iterations + kappa.\n");
            printf("for four args, add file_in and file_out.\n");
            printf("ex: computer@username: serial_main 0.1 20 file_in.jpeg file_out.jpeg\n");
            printf("\n");
        }
        if(input_state == 1){
            printf("\n");
            printf("Using default values, kappa = 0.1, iters = 10\n");
            printf("\n");
        }
        if(input_state == 2){
            printf("\n");
            printf("Error: give two, or four arguments, order: kappa, iters, innput_file, output_file. --EXITING--\n");
            printf("\n");
        }
        if(input_state == 3){
            printf("\n");
            printf("Using default filenames\n");
            printf("\n");
        }
        if(input_state == 4){
            printf("\n");
            printf("kappa must be in the range 0-1, --EXITING--\n");
            printf("\n");
        }
        if(input_state == 5){
            printf("\n");
            printf("Error: give two, or four arguments, order: kappa, iters, innput_file, output_file. --EXITING--\n");
            printf("\n");
        }
        if(input_state == 6){
            printf("\n");
            printf("set: kappa= %f, iterations =  %d, \n input= %s, output = %s\n", atof(argv[2]), atoi(argv[1]), argv[3], argv[4]);
            printf("\n");
        }
        if(input_state == 7){
            printf("\n");
            printf("kappa must be in the range 0-1, --EXITING--\n");
            printf("\n");
        }
        if(input_state == 8){
            printf("\n");
            printf("too many arguments. max = 4, --EXITING--\n");
            printf("\n");
        }
    }

    //end of input handling

    //take care of invalid input and abort
    if (error != 0) {
        if (my_rank == 0) {
            fprintf(stderr, "Error: Program terminated with error code %d\n", error);
        }
        MPI_Finalize ();
        //MPI_Abort(MPI_COMM_WORLD,0);
        exit(error);
        //return 0;
    }

    MPI_Barrier(MPI_COMM_WORLD); //check input before continue

    if (my_rank==0){  //import data
        // get 1d data
        printf("proc %d, importing data \n", my_rank);
        import_JPEG_file(input_jpeg_filename, &image_chars, &image_height, &image_width, &num_components);
        allocate_image(&whole_image, image_height, image_width);
        convert_jpeg_to_image (image_chars, & whole_image);
    }

    MPI_Bcast(&image_width,1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&image_height,1, MPI_INT, 0, MPI_COMM_WORLD);

    //divide m x n (AND COMNQUER!)
    my_m = (int)(image_height/(num_procs-1)); //height = total_height/num of processes --------------------
    my_n = (int)image_width; //width, let width be same as before
    if(my_rank == (num_procs-1)){
        my_m = image_height - ((my_rank-1)*(my_m));
    }

    //printf("proc %d, my_n = %d, my_m = %d, number of processes = %d  \n",my_rank, my_n, my_m, num_procs); //not necessary, but gives a report of every process

    //allocate memory and fill with data (for each process)
    allocate_image(&u, my_m+6, my_n);
    allocate_image(&u_bar, my_m+6, my_n);



    /* each process asks process 0 for a partitioned region */
    /* of image_chars and copy the values into u */
    /* ... */

    for(dummy = 0; dummy< iters; dummy++){//putting iterations outside loop

        /* Here we have to make a choice wheter we let each process */
        /* comunicate inbetween each iteration or compute every iteration */
        /* and then comunicate the result. I chose to go with the former */
        /* to ensure that we compute correctly. the latter will give allot less */
        /* overhead though. And therefore be the "quickest" solution */

        if(my_rank != 0){
            for(i = 0; i< my_m+6; i++){
                MPI_Recv(& u.image_data[i][0], my_n , MPI_FLOAT,0, 100,MPI_COMM_WORLD, &status);

            }
        }
        if (my_rank == 0){
            for(k = 0; k< my_m+6; k++){  //first process handeled uniquely
                MPI_Send(& whole_image.image_data[k][0], my_n , MPI_FLOAT, 1, 100, MPI_COMM_WORLD);

            }
            for(i = 1; i<= num_procs-3; i++){//everything inbetween
                for(k = (i*my_m)-3; k< ((i*my_m) + my_m) +3; k++){
                    MPI_Send(& whole_image.image_data[k][0], my_n , MPI_FLOAT, i+1, 100, MPI_COMM_WORLD);
                }
            }
            for(k = my_m*(num_procs-2)-6; k< image_height; k++){  //last process handeled uniquely
                MPI_Send(& whole_image.image_data[k][0], my_n , MPI_FLOAT, num_procs-1, 100, MPI_COMM_WORLD);
            }
        }

        /*----------------------------------------- */
        /* iso diffusion                            */
        /*----------------------------------------  */


        if(my_rank != 0){
            //smoooothen
            iso_diffusion_denoising (&u, &u_bar, kappa, 1); //putting iterations outside loop
        }


        /* each process sends its resulting content of u_bar to process 0     */
        /* process 0 receives from each process incoming values and */
        /* copy them into the designated region of struct whole_image */
        /* ... */


        if(my_rank != 0){
            if(my_rank == 1){//first partitioned area
                //printf ("proc %d sending diffused \n", my_rank);
                for(i = 0; i<= my_m+2; i++){
                    MPI_Send(& u.image_data[i][0], my_n, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
                }
            }else if(my_rank == num_procs-1){//last partitioned area
                for(i = 0; i< my_m; i++){
                    MPI_Send(& u.image_data[i+6][0], my_n, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
                }
            }else{//inbetween
                for(i = 0; i< my_m; i++){
                    MPI_Send(& u.image_data[i+3][0], my_n, MPI_FLOAT, 0, 100, MPI_COMM_WORLD);
                }
            }
        }

        if (my_rank == 0){
            for(k = 0; k<= my_m+2; k++){  //first process handeled uniquely
                MPI_Recv(& whole_image.image_data[k][0], my_n , MPI_FLOAT, 1, 100, MPI_COMM_WORLD, &status);
            }
            for(i = 1; i<= num_procs-3; i++){//everything inbetween
                for(k = (i*my_m); k< ((i*my_m) + my_m); k++){
                    MPI_Recv(& whole_image.image_data[k][0], my_n , MPI_FLOAT, i+1, 100, MPI_COMM_WORLD, &status);
                }
            }
            for(k = my_m*(num_procs-2); k< image_height; k++){  //last process handeled uniquely
                MPI_Recv(& whole_image.image_data[k][0], my_n , MPI_FLOAT, num_procs-1, 100, MPI_COMM_WORLD, &status);
            }
        }
    }

    //2d ---> 1d
    if (my_rank==0) {
        i = 0;
        printf("proc %d -- exporting image \n", my_rank);
        convert_image_to_jpeg(&whole_image, image_chars);
        export_JPEG_file(output_jpeg_filename, image_chars, image_height, image_width, num_components, 75);
    }

    MPI_Barrier (MPI_COMM_WORLD);//carefull not to free annything before used #segfault

    //free memory
    if (my_rank==0) {
        deallocate_image (&whole_image);
        printf("done");
    }

    deallocate_image(&u);
    deallocate_image(&u_bar);
    MPI_Finalize ();
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
    /* but this would give allot of overhead */

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



