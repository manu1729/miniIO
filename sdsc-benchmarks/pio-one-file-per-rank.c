/*
 * Testing POSIX I/O
 * writes one file per rank. Basically writes all array variables 
 * into a single file one after the other, avoiding significant number of 
 * file open / close operations. However, there is fwrite call per variable.
 *
*/

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<unistd.h>
#include<getopt.h>
#include "pio.h"

int main(int argc, char * argv[]) {

 int rank, size, filesize=-1;
 int darrsize, iarrsize,blk;

 char fname[100] = {'\0'};
 double *arr[16];
 
 MPI_Init(&argc, &argv);
 MPI_Comm_size(MPI_COMM_WORLD, &size);
 MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 int st = read_params(argc, argv, &filesize, fname, &blk,rank);

if (st == -1) {
   MPI_Finalize();
   return -1;
 }

 if (filesize == -1)
   filesize = FILE_SIZE_PER_VAR;
 int i,j;
 if (fname[0] == '\0') {
  sprintf(fname, "%s/test-fall.%d", FILE_PATH,rank);
 }
 for (i = 0; i < 16; i++) {
  arr[i] = (double *) malloc(filesize*1024);
 }



// initialize the arrays with random values;
 srand(time(NULL));
 for (j = 0; j < 16; j++) {
  for (i = 0; i < ((filesize*1024)/sizeof(double)); i++) {
   arr[j][i] = rand() * 1.0;
  }
 }
   
 // write out the arrays into files
 double s =  MPI_Wtime();
 FILE *fp = fopen(fname, "wb+");
 for (j=0;j<blk;j++) {
   for (i=0;i < 16; i++) {
    fwrite(arr[i],sizeof(double),((filesize*1024)/sizeof(double)),fp);
   }
 }
 for (i=0;i < 16; i++)
  free(arr[i]);
 fclose(fp);
 MPI_Barrier(MPI_COMM_WORLD);

 double e =  MPI_Wtime();
 if (!rank)
   printf("Time: %lf\n", (e-s));
 MPI_Finalize();

 return 0;
}
