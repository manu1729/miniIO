/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <mpi.h>

#include "hdf5.h"
#include "hdf5struct.h"

#define TIMEIO

static const int fnstrmax = 4095;

void writehdf5_quilt(const int num_varnames, char varnames[][20], MPI_Comm comm, int tstep,
               int ni, int nj, int nk, int cni, int cnj, int cnk, 
               float *data, float *height, int *ola_mask, int *ol_mask, hsize_t *h5_chunk, int hdf5_compress, int rank)
{
    char fname[fnstrmax+1];
    char fname_xdmf[fnstrmax+1];
    int timedigits = 4;
    MPI_Info info = MPI_INFO_NULL;

    hid_t file_id;
    hid_t plist_id;
    hid_t memspace;
    hid_t filespace;
    hid_t did;
    hsize_t start[3], count[3], stride[3];
    hsize_t *block=NULL;
    hsize_t dimsf[3];
    hsize_t dimsm[3];
    int j;
    herr_t err;
    hid_t chunk_pid;
    hsize_t chunk[3];
    
    snprintf(fname, fnstrmax, "struct_t%0*d.h5", timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "struct_t%0*d.xmf", timedigits, tstep);

    dimsf[0] = nk;
    dimsf[1] = nj;
    dimsf[2] = ni;
    dimsm[0] = cnk;
    dimsm[1] = cnj;
    dimsm[2] = cni;

#ifdef TIMEIO
    double createfile, prewrite, write, postwrite;   /* Timers */
    timer_tick(&createfile, comm, 0);
#endif

  if (rank == 0) {
      if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5 error: could not create %s \n", fname);
	MPI_Abort(comm, 1);
      }

      filespace = H5Screate_simple(3, dimsf, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	if( H5Pset_fill_time(chunk_pid, H5D_FILL_TIME_NEVER) < 0 ) {
	  printf("writehdf5 error: Could not set fill time\n");
	  MPI_Abort(comm, 1);
	}
	chunk[0]=dimsm[0]/h5_chunk[1];
	chunk[1]=dimsm[1]/h5_chunk[0];
	chunk[2]=dimsm[2];
	
	H5Pset_chunk(chunk_pid, 3, chunk);

	if(hdf5_compress == 1) {
	  
	  /* Set ZLIB / DEFLATE Compression using compression level 6. */
	  H5Pset_deflate (chunk_pid, 6);
	  
	  /* Uncomment these lines to set SZIP Compression
	     szip_options_mask = H5_SZIP_NN_OPTION_MASK;
	     szip_pixels_per_block = 16;
	     status = H5Pset_szip (plist_id, szip_options_mask, szip_pixels_per_block);
	  */
	}
      }

      /* Create the dataset with default properties */
      
      for (j=0; j<num_varnames; j++) {

	if(strcmp(varnames[j],"data") == 0 || strcmp(varnames[j],"height") == 0) {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
	} else {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_INT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
	}

	H5Dclose(did);
      }

      /* Close filespace. */
      H5Sclose(filespace);
      if(h5_chunk) {
	/* Close . */
	H5Pclose(chunk_pid);
      }
      /* Close the file */
      H5Fclose(file_id);

      /* Create xdmf file for timestep */
//      write_xdmf_xml(fname, fname_xdmf, num_varnames, varnames, ni, nj, nk, deltax, deltay, deltaz);

  }

#ifdef TIMEIO
    timer_tock(&createfile);
    timer_collectprintstats(createfile, comm, 0, "CreateFile");
    timer_tick(&prewrite, comm, 0);
#endif

    /* Set up MPI info */
//    MPI_Info_create(&info);
//    MPI_Info_set(info, "striping_factor", "1");

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    if(H5Pset_fapl_mpio(plist_id, comm, info) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_all_coll_metadata_ops(plist_id, 1);
    H5Pset_coll_metadata_write(plist_id, 1);

    if( (file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id)) < 0) {
      fprintf(stderr, "writehdf5p error: could not open %s \n", fname);
      MPI_Abort(comm, 1);
    }

    if(H5Pclose(plist_id) < 0) {
      printf("writehdf5p error: Could not close property list \n");
      MPI_Abort(comm, 1);
    }

    int is, js, ks;
    is = ola_mask[cni*cnj*cnk] * cni;
    js = ola_mask[cni*cnj*cnk+1] * cnj;
    ks = ola_mask[cni*cnj*cnk+2] * cnk;

    if(h5_chunk) {
      block = malloc(3 * sizeof(hsize_t));
      count[0] = 1;
      count[1] = dimsm[1];
      count[2] = dimsm[2];
      start[0] = (hsize_t)(ks);
      start[1] = (hsize_t)(js);
      start[2] = (hsize_t)(is);
      block[0] = dimsm[0];
      block[1] = 1;
      block[2] = 1;
    } else {
      start[0] = (hsize_t)(ks);
      start[1] = (hsize_t)(js);
      start[2] = (hsize_t)(is);
      count[0] = (hsize_t)(cnk);
      count[1] = (hsize_t)(cnj);
      count[2] = (hsize_t)(cni);
    }

    memspace = H5Screate_simple(3, dimsm, NULL);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

#ifdef TIMEIO
    timer_tock(&prewrite);
    timer_collectprintstats(prewrite, comm, 0, "PreWrite");
    
#endif
    for (j=0; j<num_varnames; j++) {
#ifdef TIMEIO
     timer_tick(&write, comm, 0);
#endif
      did = H5Dopen(file_id, varnames[j],H5P_DEFAULT);

      /* 
       * Each process defines dataset in memory and writes it to the hyperslab
       * in the file.
       */

      /* Select hyperslab in the file.*/
      filespace = H5Dget_space(did);
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, block);
      
      err = 0;
      if(strcmp(varnames[j],"data") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, data);
      } else if(strcmp(varnames[j],"height") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace , plist_id, height);
      } else if(strcmp(varnames[j],"ola_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, ola_mask);
      } else if(strcmp(varnames[j],"ol_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, ol_mask);
      } else {
	printf("writehdf5 error: Unknown how to handle variable %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }

      if( err < 0) {
	fprintf(stderr, "writehdf5 error: could not write datset %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }
      err = H5Dclose(did);

#ifdef TIMEIO
    timer_tock(&write);
    timer_collectprintstats(write, comm, 0, varnames[j]);
#endif
    }
#ifdef TIMEIO
    timer_tick(&postwrite, comm, 0);
#endif
    if(H5Sclose(memspace)  != 0)
      printf("writehdf5 error: Could not close memory space \n");

    if(H5Pclose(plist_id) < 0)
      printf("writehdf5 error: Could not close property list \n");

    if(h5_chunk) {
      free(block);
    }

    if(H5Fclose(file_id) != 0)
      printf("writehdf5 error: Could not close HDF5 file \n");
#ifdef TIMEIO
    timer_tock(&postwrite);
    timer_collectprintstats(postwrite, comm, 0, "PostWrite");
#endif
}

void test(float *data, ...) {
  printf("inside test ****************\n");
}

void writehdf5_quilt_new(const int num_varnames, char varnames[][20], MPI_Comm comm, int tstep,
               int ni, int nj, int nk, int cni, int cnj, int cnk, 
               float *data, float *height, int *ola_mask, int *ol_mask, hsize_t *h5_chunk, int hdf5_compress, int rank)
{
    char fname[fnstrmax+1];
    char fname_xdmf[fnstrmax+1];
    int timedigits = 4;
    MPI_Info info = MPI_INFO_NULL;

    hid_t file_id;
    hid_t plist_id;
    hid_t memspace;
    hid_t filespace;
    hid_t did;
    hsize_t start[3], count[3], stride[3];
    hsize_t *block=NULL;
    hsize_t dimsf[3];
    hsize_t dimsm[3];
    int j;
    herr_t err;
    hid_t chunk_pid;
    hsize_t chunk[3];
    
    snprintf(fname, fnstrmax, "struct_t%0*d.h5", timedigits, tstep);
    snprintf(fname_xdmf, fnstrmax, "struct_t%0*d.xmf", timedigits, tstep);

    dimsf[0] = nk;
    dimsf[1] = nj;
    dimsf[2] = ni;
    dimsm[0] = nk;
    dimsm[1] = nj;
    dimsm[2] = cni;

#ifdef TIMEIO
    double createfile, prewrite, write, postwrite;   /* Timers */
    timer_tick(&createfile, comm, 0);
#endif
  int r;
  MPI_Comm_rank(comm, &r);
  if (!r) {
      if( (file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0) {
	fprintf(stderr, "writehdf5 error: could not create %s \n", fname);
	MPI_Abort(comm, 1);
      }

      filespace = H5Screate_simple(3, dimsf, NULL);

      chunk_pid = H5Pcreate(H5P_DATASET_CREATE);
      if(h5_chunk) {
	H5Pset_layout(chunk_pid, H5D_CHUNKED);
	if( H5Pset_fill_time(chunk_pid, H5D_FILL_TIME_NEVER) < 0 ) {
	  printf("writehdf5 error: Could not set fill time\n");
	  MPI_Abort(comm, 1);
	}
	chunk[0]=dimsm[0]/h5_chunk[1];
	chunk[1]=dimsm[1]/h5_chunk[0];
	chunk[2]=dimsm[2];
	
	H5Pset_chunk(chunk_pid, 3, chunk);

	if(hdf5_compress == 1) {
	  
	  /* Set ZLIB / DEFLATE Compression using compression level 6. */
	  H5Pset_deflate (chunk_pid, 6);
	  
	}
      }

      /* Create the dataset with default properties */
      
      for (j=0; j<num_varnames; j++) {

	if(strcmp(varnames[j],"data") == 0 || strcmp(varnames[j],"height") == 0) {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_FLOAT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
	} else {
	  did = H5Dcreate(file_id, varnames[j], H5T_NATIVE_INT, filespace, H5P_DEFAULT, chunk_pid, H5P_DEFAULT);
	}

	H5Dclose(did);
      }

      /* Close filespace. */
      H5Sclose(filespace);
      if(h5_chunk) {
	/* Close . */
	H5Pclose(chunk_pid);
      }
      /* Close the file */
      H5Fclose(file_id);

  }
  MPI_Barrier(comm);
#ifdef TIMEIO
    timer_tock(&createfile);
    timer_collectprintstats(createfile, comm, 0, "CreateFile");
    timer_tick(&prewrite, comm, 0);
#endif

    /* Set up MPI info */
    MPI_Info_create(&info);
    MPI_Info_set(info, "striping_factor", "1");

    /* Set up file access property list with parallel I/O access */
    if( (plist_id = H5Pcreate(H5P_FILE_ACCESS)) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_libver_bounds(plist_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);

    if(H5Pset_fapl_mpio(plist_id, comm, info) < 0) {
      printf("writehdf5p error: Could not create property list \n");
      MPI_Abort(comm, 1);
    }

    H5Pset_all_coll_metadata_ops(plist_id, 1);
    H5Pset_coll_metadata_write(plist_id, 1);

    if( (file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id)) < 0) {
      fprintf(stderr, "writehdf5p error: could not open %s \n", fname);
      MPI_Abort(comm, 1);
    }

    if(H5Pclose(plist_id) < 0) {
      printf("writehdf5p error: Could not close property list \n");
      MPI_Abort(comm, 1);
    }

    int is, js, ks;

    memspace = H5Screate_simple(3, dimsm, NULL);

    /* Create property list for collective dataset write. */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

#ifdef TIMEIO
    timer_tock(&prewrite);
    timer_collectprintstats(prewrite, comm, 0, "PreWrite");
    
#endif
    int i;
    int *t = ola_mask;
    float *xfordataf = (float *) malloc(cni*nj*nk*sizeof(float));
    int *xfordatai = (int *) malloc(cni*nj*nk*sizeof(int));
    for (j=0; j< num_varnames; j++) {

#ifdef TIMEIO
     timer_tick(&write, comm, 0);
#endif
      did = H5Dopen(file_id, varnames[j],H5P_DEFAULT);
      filespace = H5Dget_space(did);
      printf("%d ---------------%s\n", r, varnames[j]);
      if(strcmp(varnames[j],"data") == 0)
        xform_dataf(r,xfordataf,data,ni,nj,nk,cni,cnj,cnk,ola_mask);
      else if(strcmp(varnames[j],"height") == 0) 
        xform_dataf(r,xfordataf,height,ni,nj,nk,cni,cnj,cnk,ola_mask);
      else if(strcmp(varnames[j],"ol_mask") == 0) 
        xform_datai(r,xfordatai,ol_mask,ni,nj,nk,cni,cnj,cnk,ola_mask,0);
      else if(strcmp(varnames[j],"ola_mask") == 0) 
        xform_datai(r,xfordatai,ola_mask,ni,nj,nk,cni,cnj,cnk,ola_mask,1); 

      start[0] = 0;
      start[1] = 0;
      start[2] = cni*r;
      count[0] = nk;
      count[1] = nj;
      count[2] = cni;
    
  
      /* Select hyperslab in the file.*/
      H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, block);
        
      err = 0;
      if(strcmp(varnames[j],"data") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace, plist_id, xfordataf);
      } else if(strcmp(varnames[j],"height") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_FLOAT, memspace, filespace , plist_id, xfordataf);
      } else if(strcmp(varnames[j],"ola_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, xfordatai);
      } else if(strcmp(varnames[j],"ol_mask") == 0) {
	err = H5Dwrite(did, H5T_NATIVE_INT, memspace, filespace, plist_id, xfordatai);
      } else {
	printf("writehdf5 error: Unknown how to handle variable %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }

      if( err < 0) {
	fprintf(stderr, "writehdf5 error: could not write datset %s \n", varnames[j]);
	MPI_Abort(comm, 1);
      }
      err = H5Dclose(did);

#ifdef TIMEIO
    timer_tock(&write);
    timer_collectprintstats(write, comm, 0, varnames[j]);
#endif
    }
#ifdef TIMEIO
    timer_tick(&postwrite, comm, 0);
#endif
    if(H5Sclose(memspace)  != 0)
      printf("writehdf5 error: Could not close memory space \n");

    if(H5Pclose(plist_id) < 0)
      printf("writehdf5 error: Could not close property list \n");

    if(h5_chunk) {
      free(block);
    }

    if(H5Fclose(file_id) != 0)
      printf("writehdf5 error: Could not close HDF5 file \n");
#ifdef TIMEIO
    timer_tock(&postwrite);
    timer_collectprintstats(postwrite, comm, 0, "PostWrite");
#endif
}


void xform_dataf(int rank, float *xformdata, float *data, int ni, int nj,int nk, int cni,int cnj,int cnk, int *offs) {
  int l,i,j,k,is,js,ks;
  int offset, *t, c=0;
  for (l = 0; l < nj*nk/(cnj*cnk)  /*temporary fix, should be #compute-tasks*/; l++) {
    offset = l*cni*cnj;
    for (k=0;k<cnk;k++)
     for (j=0;j<cnj;j++)
      for (i=0;i<cni;i++) {
        xformdata[offset + cni*nj*k + cni*j + i] = data[c++];
//        printf("%d %f \n",offset + cni*nj*k + cni*j + i, data[c-1]);
      }
  }
}

void xform_datai(int rank, int *xformdata, int *data, int ni, int nj,int nk, int cni,int cnj,int cnk, int *offs, int flag) {
  int l,i,j,k,is,js,ks;
  int offset, *t, c=0;
  for (l = 0; l < nj*nk/(cnj*cnk) /*temporary fix, should be #compute-tasks*/; l++) {
    offset = l*cni*cnj;
    for (k=0;k<cnk;k++)
     for (j=0;j<cnj;j++)
      for (i=0;i<cni;i++)
        xformdata[offset + cni*nj*k + cni*j + i] = data[c++];
    if (flag) c += 3;
  }
}

