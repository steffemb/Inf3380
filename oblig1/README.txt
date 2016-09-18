------------------------------------------------------------
USAGE:

Both programs can take 0, 2, or 4 arguments. If two are given they
should be "iterations" and "kappa". If four are given they
should be "iterations", "kappa", "infile", and "outfile" 

serial example: 
$ ./serial_main 20 0.2 infile.jpeg outfile.jpeg

paralell example: 
$ mpirun -n 4  ./parallel_main 20 0.2 infile.jpeg outfile.jpeg

If none or two arguments are given, the program will use 
default values.


------------------------------------------------------------



To complete oblig 1 you have to first write a serial program
by filling out the serial_main.c file in the serial/
directory and finally write a parallel program by filling
out the parallel_main.c file in the parallel/ directory.

The two programs shall perform an iso diffusion denoising on
a JPEG image. The two programs shall also accept the
following parameters (in order)

$ ./program number_of_iterations kappa_value infile outfile

When you are in the serial/ or parallel/ directories you may
compile the programs by typing "make" in the terminal.

When you are ready to deliver you enter the base directory
and type "make delivery" in the terminal. A tarball with
your updated files will then be created ready for delivery.
