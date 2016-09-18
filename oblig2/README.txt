DESCRIPTION:

This is a program that is ment to answer the task in "obligatory assignment 2"
of the course "inf3380 vår2016".

the program uses a column-vise block partition (yes i know we should have used 
row-wise, but i asked during group session and got a thumbs up) on the 
MPI processes and has OpenMP implemented inside the function 
matrix_multiplication wich each process calls.

USAGE:

 program takes 0,1 or two arguments. 
 0 arguments --> default values

 1 argument --> command

 list of valid commands: 

 -print_matrix
 -bench
 -check_error

 the last option is to specify input and output filenames:
 example: user@pc~/destination_folder/$ mpirun -np 4 ./matmultiply_main -bench

ADDED FILES:

-a MAKE file is added for -make and -make clean (only tested on 32-bit linux).
-the file matrix_c.bin is a copy of small_matrix_c.bin from the course site.
 It is used for reference checking with the command -check_error.
-The rest of the filesare also from the course site exept my_matrix_c.bin
 this is the written file from the program. 


WARNING: number of MPI processes must be dividable by 2!
