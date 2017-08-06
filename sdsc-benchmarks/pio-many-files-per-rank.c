/*
 * Testing POSIX I/O
 * writes one file per each of the 16 array variables and per rank,
 * so there will be plenty of small files being written
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
 int darrsize, iarrsize,blk=1;

 char *fname[16] = {NULL};
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
 if (fname[0] == NULL) {
  for (i = 0; i < 16; i++) {
   fname[i] = (char*) malloc(100);
   sprintf(fname[i], "%s/test-f%d.%d", FILE_PATH,i,rank);
  }
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
 FILE *fp;
 for (j=0; j< blk; j++) {
   for (i=0;i < 16; i++) {
    if (j == 0)
      fp = fopen(fname[i], "wb+");
    else 
      fp = fopen(fname[i], "ab+");
    fwrite(arr[i],sizeof(double),((filesize*1024)/sizeof(double)),fp);
    fclose(fp);
   }
 }
 for (i=0;i < 16; i++)
   free(arr[i]);

 MPI_Barrier(MPI_COMM_WORLD);

 double e =  MPI_Wtime();
 if (!rank)
   printf("Time: %lf\n", (e-s));
 MPI_Finalize();

 return 0;
}
