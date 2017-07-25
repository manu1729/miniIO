/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/
#ifdef HAS_HDF5
#  include "hdf5.h"
#endif
void writehdf5(const int num_varnames, char **varnames, MPI_Comm comm, int rank, int nprocs, int tstep, 
	       int is, int js, int ks,
               int ni, int nj, int nk, int cni, int cnj, int cnk, 
               float deltax, float deltay, float deltaz, 
               float *data, float *height, int *ola_mask, int *ol_mask, hsize_t *h5_chunk, int hdf5_compress);
void writehdf5_quilt(const int num_varnames, char varnames[][20], MPI_Comm comm, int tstep,
               int ni, int nj, int nk, int cni, int cnj, int cnk,
               float *data, float *height, int *ola_mask, int *ol_mask, hsize_t *h5_chunk, int hdf5_compress,int i);

void writehdf5_quilt_test(const int num_varnames, char varnames[][20], MPI_Comm comm, int tstep,
               int ni, int nj, int nk, int cni, int cnj, int cnk,
               float *data, float *height, int *ola_mask, int *ol_mask, hsize_t *h5_chunk, int hdf5_compress,int i);

void test(float *, ...);
void xform_dataf(int rank, float *xformdata, float *data, int ni, int nj,int nk, int cni,int cnj,int cnk, int *offs);
void xform_datai(int rank, int *xformdata, int *data, int ni, int nj,int nk, int cni,int cnj,int cnk, int *offs, int flag);
