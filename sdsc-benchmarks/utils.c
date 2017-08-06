/*
 * Utility functions
 * 
 *
*/

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<getopt.h>
#include "utils.h"

int read_params(int argc, char *argv[], int *fs, char *fname, int *blk, int rank) {

 int opt = 0;
 while (( opt = getopt(argc, argv, "s:p:b:")) != -1) {
   switch (opt) {
     case 'b': *blk = atoi(optarg);
       break;
     case 's': *fs = atoi(optarg);
       break;
     case 'p': sprintf(fname, "%s/test-fall%d", optarg,rank);
       break;
     default:  printf("Error is arguments\n");
       return -1;
   }
 }
 return 0;

}
