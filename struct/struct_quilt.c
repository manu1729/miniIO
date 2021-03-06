/*
* Copyright (c) DoD HPCMP PETTT.  All rights reserved.
* See LICENSE file for details.
*/

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>
#include <mpi.h>
#include "open-simplex-noise.h"
#include "timer.h"


/* #include <limits.h> */
/* #include <assert.h> */

static const float FILLVALUE = -999;

#ifdef HAS_ADIOS
#  include "adiosstruct.h"
#endif

#ifdef HAS_HDF5
#  include "hdf5struct.h"
#endif

void print_usage(int rank, const char *errstr);

int main(int argc, char **argv)
{
  int debug=0;                  /* Flag to generate debug prints statements */
  int debugIO=0;                /* Flag to generate debug IO*/  
  int a, i, j, k, t;            /* loop indices */
  int tt;                       /* Actual time step from tstart */
  float x, y, z;
  double noisespacefreq = 10.0; /* Spatial frequency of noise */
  double noisetimefreq = 0.25;  /* Temporal frequency of noise */
  int tstart = 0;
  int nt = 10;                  /* Number of time steps */
  int ni = 0;                   /* Global grid size */
  int nj = 0;
  int nk = 0;
  int inp = 0;      /* Number of tasks in i */
  int jnp = 0;      /* Number of tasks in j */
  int knp = 1;      /* Number of tasks in k */
  int iots = 0;     /* Number of IO tasks: Manu*/
  int numtask;                   /* task per processor */
  int point_id, tmp_id;          /* grid point id */
  int x_index, y_index, z_index; /* point index along each axis */
  int xy_dims,  x_dims;
  float deltax, deltay, deltaz; 
  int numPoints;
  float *data;
  float *rowdata; /* Added for quilting - Manu */
  float *height;
  float *rowheight; /* Added for quilting - Manu */
  int height_index;
  int hindex;
  int maskTindex;
  int *ola_mask;
  int *ol_mask;
  int *rowola_mask; /* Added for quilting - Manu */
  int *rowol_mask; /* Added for quilting - Manu */
  float mask_thres=0.0;     /* upper mask threshold  (range -1 to 1) */
  float bot_mask_thres=0.5; /* bottom mask threshold (range 0.0 to (mask_thres+1)/2 ) */
  int mask_thres_index;
  struct osn_context *simpnoise;    /* Open simplex noise context */
  double heighttime, computetime, outtime;   /* Timers */
  double apptime;   /* Timers */
  
  const int num_varnames=4;
  char varnames[4][20] = {"data", "height", "ola_mask", "ol_mask"};

  /* MPI vars */
  MPI_Comm comm = MPI_COMM_WORLD;
  int cprocs[3], cpers[3], crnk[3] = {-1,-1,-1};  /* MPI Cartesian info */
  int rank, nprocs;
  int cni, cnj, cnk;   /* Points in this task */
  int is, js, ks;     /* Global index starting points */
  float xs, ys, zs;    /* Global coordinate starting points */  

  /* MPI vars for IO quilting: Manu */
  int new_rank;
  MPI_Comm  new_comm;
  MPI_Comm inter_comm, row_comm;
  MPI_Comm inter_comm_rows[20];
  int color;
 
  /* ADIOS vars */
  uint64_t cstart=0;
  uint64_t cnpoints=0;
  uint64_t npoints=0;

#ifdef HAS_HDF5
  int hdf5out = 0;
  hsize_t *hdf5_chunk=NULL;
  int hdf5_compress = 0;
#endif

#ifdef HAS_ADIOS
  char      *adios_groupname="struct";
  char      *adios_method=NULL;   /* POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5   */
  struct adiosstructinfo adiosstruct_nfo;
#endif

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  /* Parse command line */
  for(a = 1; a < argc; a++) {

    if(!strcasecmp(argv[a], "--tasks")) {
            inp = atoi(argv[++a]);
            jnp = atoi(argv[++a]);
    } else if (!strcasecmp(argv[a], "--iotasks")) {   /* For quilting: Manu*/
            iots = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--size")) {
      ni = atoi(argv[++a]);
      nj = atoi(argv[++a]);
      nk = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--maskthreshold")) {
      mask_thres = strtof(argv[++a], NULL);
      bot_mask_thres = (mask_thres+1.0)/2;
    } else if(!strcasecmp(argv[a], "--noisespacefreq")) {
      noisespacefreq = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--noisetimefreq")) {
      noisetimefreq = strtod(argv[++a], NULL);
    } else if(!strcasecmp(argv[a], "--tsteps")) {
      nt = atoi(argv[++a]);
    } else if(!strcasecmp(argv[a], "--tstart")) {
      tstart = atoi(argv[++a]);
    }else if(!strcasecmp(argv[a], "--debug")) {
      debug = 1;
    }
#ifdef HAS_ADIOS
    else if(!strcasecmp(argv[a], "--adios")) {
      adios_method = argv[++a];
    }
    else if(!strcasecmp(argv[a], "--debugIO")) {
      debugIO = 1;
    }
#endif
    else if( !strcasecmp(argv[a], "--hdf5") )
      {
#ifdef HAS_HDF5
	hdf5out = 1;
      }
    else if(!strcasecmp(argv[a], "--hdf5_chunk")) {
      hdf5_chunk = malloc(2 * sizeof(hsize_t));
      hdf5_chunk[1] = (hsize_t)strtoul(argv[++a], NULL, 0);
      hdf5_chunk[0] = (hsize_t)strtoul(argv[++a], NULL, 0);
      if (hdf5_chunk[0] <= 0 || hdf5_chunk[1] <= 0 ) {
	print_usage(rank, "Error: Illegal chunk dim sizes");
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    else if(!strcasecmp(argv[a], "--hdf5_compress")) {
      hdf5_compress=1;
    }
#else
      if(rank == 0)   fprintf(stderr, "HDF5 option not available: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
#endif
    else {
      if(rank == 0)   fprintf(stderr, "Option not recognized: %s\n\n", argv[a]);
      print_usage(rank, NULL);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
   }

  numPoints = ni*nj*nk;
  npoints  = numPoints;
 
  /* Check arguments & proc counts */
  if(inp < 1 || jnp < 1 ) {
    print_usage(rank, "Error: tasks not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(ni < 1 || nj < 1 || nk < 1) {
    print_usage(rank, "Error: size not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if((inp*jnp + iots) != nprocs) {   /* For quilting: Manu */
    print_usage(rank, "Error: number of required tasks does not equal total MPI tasks");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  if(inp != iots) {   /* For quilting: Manu */
    print_usage(rank, "Error: number of I/O tasks should be equal to the number of compute tasks along i-drection");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(ni % inp || nj % jnp || nk % knp) {
    print_usage(rank, "Error: number of points on an axis is not evenly divisible "
                "by axis tasks.\n   This is required for proper load balancing.");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if(nt < 1 ) {
    print_usage(rank, "Error: number of timesteps not specified or incorrect");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (iots == 0) {
     /* Set up Cartesian communicator */
    cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
    cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */
    MPI_Cart_create(MPI_COMM_WORLD, 3, cprocs, cpers, 1, &comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Cart_coords(comm, rank, 3, crnk);
  } else {
    color = rank/(inp*jnp*knp);
    MPI_Comm_split(MPI_COMM_WORLD,color, rank, &new_comm);
    if (rank < inp*jnp*knp) {
     cprocs[0] = inp;  cprocs[1] = jnp;  cprocs[2] = knp;
     cpers[0] = 0;  cpers[1] = 0;  cpers[2] = 0;    /* No periodicity */
     MPI_Cart_create(new_comm, 3, cprocs, cpers, 0, &comm);
     int old_rank = rank;
     MPI_Comm_rank(comm, &rank);
     MPI_Cart_coords(comm, rank, 3, crnk);
     color = crnk[0];
     MPI_Comm_split(new_comm, color, rank, &row_comm);
    }

    if (crnk[0] != -1) {
      MPI_Intercomm_create(row_comm, 0, MPI_COMM_WORLD, inp*jnp*knp+crnk[0], crnk[0], &inter_comm_rows[crnk[0]]);
    } else if (crnk[0] == -1) {
      int nr;
      MPI_Comm_rank(new_comm, &nr);
      MPI_Intercomm_create(new_comm, nr, MPI_COMM_WORLD, nr*jnp*knp, nr, &inter_comm_rows[nr]);
    }
  }

    deltax = 1.f/(ni-1);
    deltay = 1.f/(nj-1);
    deltaz = 1.f/(nk-1);
    cni = ni / inp;
    cnj = nj / jnp;
    cnk = nk / knp;
    is = crnk[0] * cni;
    js = crnk[1] * cnj;
    ks = crnk[2] * cnk;
    xs = is * deltax;
    ys = js * deltay;
    zs = ks * deltaz;
  
    xy_dims = ni * nj;
    x_dims = ni;


    /* adjust mask threshold  to compensate by bottom threshold */
    mask_thres = mask_thres - bot_mask_thres;
    mask_thres_index = (int) ( ((mask_thres+1)/2) * (nk-1));
    maskTindex = nk-1;
    
    if (debug) {
      /* printf("rank_cord(%d,%d,%d) rank=%d: %d of %d\n", crnk[0], crnk[1], crnk[2] , rank, rank+1, nprocs); */
      printf("Grid size= (%d x %d x %d) = %llu points\n", ni, nj, nk, npoints);
      /* printf("Grid deltas= (%f x %f x %f)\n", deltax, deltay, deltaz); */
      printf("Mask theshold = %f   mask threshold index = %d\n", mask_thres, mask_thres_index);
    }
    /* Set up osn */
    open_simplex_noise(12345, &simpnoise);   /* Fixed seed, for now */
  
    /* Allocate arrays */
   if (rank < inp*jnp*knp) {  // Manu
    data = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
    height = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float));
    ola_mask = (int *) malloc((size_t)(cni*cnj*cnk + 3)*sizeof(int));
    ol_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int));
    ola_mask[cni*cnj*cnk] = crnk[0]; ola_mask[cni*cnj*cnk+1] = crnk[1];ola_mask[cni*cnj*cnk+2] = crnk[2];
   } else {
    rowdata = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float)*inp);
    rowheight = (float *) malloc((size_t)cni*cnj*cnk*sizeof(float)*inp);
    rowola_mask = (int *) malloc((size_t)(cni*cnj*cnk+3)*sizeof(int)*inp);
    rowol_mask = (int *) malloc((size_t)cni*cnj*cnk*sizeof(int)*inp);
   }

    /* init ADIOS */
  #ifdef HAS_ADIOS
  
    if (adios_method) {
      adiosstruct_init(&adiosstruct_nfo, adios_method, adios_groupname, comm, rank, nprocs, nt,
  		     ni, nj, nk, is, cni, js, cnj, ks, cnk, deltax, deltay, deltaz, FILLVALUE);
      adiosstruct_addrealxvar(&adiosstruct_nfo, varnames[0], data);
  
      if (debugIO) {
        adiosstruct_addrealxvar(&adiosstruct_nfo, varnames[1], height);
        adiosstruct_addintxvar(&adiosstruct_nfo, varnames[2], ola_mask);
        adiosstruct_addintxvar(&adiosstruct_nfo, varnames[3], ol_mask);
      }
    }
  #endif
  
    /* generate masked grid */
    /* Spatial loops */
    size_t ii;     /* data index */

  timer_tick(&apptime, MPI_COMM_WORLD, 1);

  double s = MPI_Wtime();

  if (rank < inp*jnp*knp) {
    timer_tick(&heighttime, comm, 1);
    z = zs;
    for(k = 0, ii = 0; k < cnk; k++) {
      y = ys;
      for(j = 0; j < cnj; j++) {
        x = xs;
        for(i = 0; i < cni; i++, ii++) {
  	x_index = (int) (x/deltax);
  	y_index = (int) (y/deltay);
  	z_index = (int) (z/deltaz);
  
  	/* calculate point index */
  	point_id = (z_index * xy_dims) + (y_index * x_dims) + x_index;
  
  	/* Get height and subtract bottom threshold */
  	height[ii] =  (float)open_simplex_noise2(simpnoise, x*noisespacefreq, y*noisespacefreq)  - bot_mask_thres;
  
  	/* height_index = (int) height[ii]/deltaz; */
  	height_index = (int) (((height[ii]+1)/2) * (nk-1));
  	  
  	/* Calculate ola_mask values */
  	if (z_index > mask_thres_index  && z_index > height_index) {
  	  ola_mask[ii] = 2;  /* Atmosphere */
  	}
  	else if (z_index < mask_thres_index  && z_index > height_index) {
  	  ola_mask[ii] = 0;  /* ocean */
  	}
  	else if (z_index <= height_index) {
  	  if (height[ii] >= mask_thres  || z_index < height_index) {
  	    ola_mask[ii] = 1;  /* land */
  	  }
  	  else {
  	    ola_mask[ii] = 0;  /* ocean */
  	  }
  	}
  	else if (z_index == mask_thres_index  && height[ii] <= mask_thres) {
  	  ola_mask[ii] = 0;  /* ocean */
  	}
  	else {
  	  printf("WARNING: ola_mask condition not considered for Point_index: (%d,%d,%d)\n"
  		 "Point_id: %d  Height: %f HeightID: %d  mask_thres_index=%d\n",
  		 x_index, y_index, z_index, point_id+1, height[ii], height_index, mask_thres_index); 
  	}
  
  	
  	hindex = (int) ( (((height[ii]+1)/2) * (nk-1)) / (mask_thres_index+1)) * (nk-1);
  
  	hindex = (int) ((height[ii]+1) / (mask_thres+1) * (nk-1));
  	
  	if (hindex > maskTindex) {
  	  hindex = maskTindex;
  	}
  
  	/* Calculate ol_mask values */
  	if (z_index < maskTindex  && z_index > hindex) {
  	  ol_mask[ii] = 0;  /* ocean */
  	}
  	else if (z_index <= hindex) {
  	  if (height[ii] >= mask_thres  || z_index < hindex) {
  	    ol_mask[ii] = 1;  /* land */
  	  }
  	  else {
  	    ol_mask[ii] = 0;  /* ocean */
  	  }
  	}
  	else if (z_index == maskTindex  && height[ii] <= mask_thres) {
  	  ol_mask[ii] = 0;  /* ocean */
  	}
  	else {
  

	  printf("WARNING: ol_mask condition not considered for Point_index: (%d,%d,%d)\n"
  		 "Point_id: %d  Height: %f HeightID: %d maskTindex=%d\n",
  		 x_index, y_index, z_index, point_id+1, height[ii], hindex, maskTindex); 
  	}
  
  	if (debug) {
  	  printf("++++++++++++++++++++++++++++++++++++++++++++\n");
  	  printf("rank_cord(%d,%d,%d) rank=%d: %d of %d\n", crnk[0], crnk[1], crnk[2] , rank, rank+1, nprocs);
  	  printf("LDims: (%d,%d,%d)\n", cni, cnj, cnk); 
  	  printf("GDims: (%d,%d,%d)\n", ni, nj, nk);
  	  printf("SDims: (%d,%d,%d)\n", is, js, ks); 
  	  printf("Point_index: (%d,%d,%d), Point_id:  %d\n", x_index, y_index, z_index,  point_id+1);
  	  printf("Point_pos: (%f, %f, %f)  mask_thres_index %d -> %d\n", x, y, z, mask_thres_index, maskTindex);
  	  printf("Height: %f HeightID: %d -> %d ola_mask=%d -> %d\n", height[ii], height_index,  hindex, ola_mask[ii], ol_mask[ii]);
  	}	      
  	
  	x += deltax;
        }
        y += deltay;
      }
      z += deltaz;
    }	
    timer_tock(&heighttime);
   }  // end of "if (rank < inp*jnp*knp)"
   
    /* generate ocean land data */
    for(t = 0, tt = tstart; t < nt; t++, tt++) {
      /* Spatial loops */
     if (rank < inp*jnp*knp) {
      timer_tick(&computetime, comm, 1);
      z = zs;
      for(k = 0, ii = 0; k < cnk; k++) {
        y = ys;
        for(j = 0; j < cnj; j++) {
  	x = xs;
  	for(i = 0; i < cni; i++, ii++) {
  
  	  if ( ol_mask[ii] == 0) {
  	    /* if  ( ola_mask[ii] == 0) { */
  	    data[ii] = (float)open_simplex_noise4(simpnoise, x*noisespacefreq, y*noisespacefreq, z*noisespacefreq, tt*noisetimefreq);
  	  }
  	  else {
  	    data[ii] = FILLVALUE;
  	  }
         
  	  x += deltax;
  	}
  	y += deltay;
        }
        z += deltaz;
      }	
      timer_tock(&computetime);
  
  
      //timer_tick(&outtime, comm, 1);
     } // end of "if (rank < inp*jnp*knp)"

    if (iots == 0)
      timer_tick(&outtime, comm, 1);
    else
      timer_tick(&outtime, new_comm, 1);
      
  #ifdef HAS_ADIOS
      if (adios_method) adiosstruct_write(&adiosstruct_nfo, tt);
  #endif
   
  #ifdef HAS_HDF5
      if(hdf5out) {
        if(rank == 0) {
  	printf("      Writing hdf5...\n");   fflush(stdout);
        }

        if (iots > 0) {
         int i;
         for (i = 0; i < inp; i++) {
           if (rank == inp*jnp*knp+i) {
            MPI_Gather(data, cni*cnj*cnk, MPI_FLOAT, rowdata, cni*cnj*cnk, MPI_FLOAT, MPI_ROOT, inter_comm_rows[i]);
           }
	 } 
         if (crnk[0] != -1)
          MPI_Gather(data, cni*cnj*cnk, MPI_FLOAT, rowdata, cni*cnj*cnk, MPI_FLOAT, crnk[0], inter_comm_rows[crnk[0]]);

         for (i = 0; i < inp; i++) {
           if (rank == inp*jnp*knp+i) {
            MPI_Gather(height, cni*cnj*cnk, MPI_FLOAT, rowheight, cni*cnj*cnk, MPI_FLOAT, MPI_ROOT, inter_comm_rows[i]);
           }
         }
         if (crnk[0] != -1)
          MPI_Gather(height, cni*cnj*cnk, MPI_FLOAT, rowheight, cni*cnj*cnk, MPI_FLOAT, crnk[0], inter_comm_rows[crnk[0]]);

         for (i = 0; i < inp; i++) {
           if (rank == inp*jnp*knp+i) {
            MPI_Gather(ola_mask, cni*cnj*cnk+3, MPI_INT, rowola_mask, cni*cnj*cnk+3, MPI_INT, MPI_ROOT, inter_comm_rows[i]);
           }
         }
         if (crnk[0] != -1)
          MPI_Gather(ola_mask, cni*cnj*cnk+3, MPI_INT, rowola_mask, cni*cnj*cnk+3, MPI_INT, crnk[0], inter_comm_rows[crnk[0]]);

         for (i = 0; i < inp; i++) {
           if (rank == inp*jnp*knp+i) {
            MPI_Gather(ol_mask, cni*cnj*cnk, MPI_INT, rowol_mask, cni*cnj*cnk, MPI_INT, MPI_ROOT, inter_comm_rows[i]);
           }
         }
         if (crnk[0] != -1)
          MPI_Gather(ol_mask, cni*cnj*cnk, MPI_INT, rowol_mask, cni*cnj*cnk, MPI_INT, crnk[0], inter_comm_rows[crnk[0]]);

        }
        if (rank >= inp*jnp*knp) {
          writehdf5_quilt_new(num_varnames, varnames, new_comm, tt, ni, nj, nk, cni, cnj, cnk, rowdata, rowheight, rowola_mask, rowol_mask, hdf5_chunk, hdf5_compress,rank);
	}

      }
  #endif
  
    timer_tock(&outtime);
    if (rank >= inp*jnp*knp) {
     timer_collectprintstats(outtime, new_comm, 0, "   Output");
    }
    if (rank < inp*jnp*knp) {  
    //  timer_tock(&outtime);
      timer_collectprintstats(computetime, comm, 0, "   Compute");
      timer_collectprintstats(outtime, comm, 0, "   Output");
      timer_collectprintstats(heighttime, comm, 0, "   Height");
    }
   }

   timer_tock(&apptime);
   timer_collectprintstats(apptime, MPI_COMM_WORLD, 0, "   Application"); 

  MPI_Barrier(MPI_COMM_WORLD);
  double e = MPI_Wtime();

  if (!rank)
   printf("Application time == %lf\n", e-s);
  
      /* finalize ADIOS */
  #ifdef HAS_ADIOS
    if (adios_method) adiosstruct_finalize(&adiosstruct_nfo);
  #endif
  
  #ifdef HAS_HDF5
    if(hdf5_chunk)
      free(hdf5_chunk);
  #endif
  
    open_simplex_noise_free(simpnoise);
  
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank < inp*jnp*knp) {
    free(data);
    free(height);
    free(ola_mask);
    free(ol_mask);
   } else {
     free(rowdata);
     free(rowheight);
     free(rowola_mask);
     free(rowol_mask); 
   }
  MPI_Barrier(MPI_COMM_WORLD);
//  MPI_Comm_free(&inter_comm);
  MPI_Finalize();

  return 0;
}


void print_usage(int rank, const char *errstr)
{
  if(rank != 0)  return;
  if(errstr)
    fprintf(stderr, "%s\n\n", errstr);
  fprintf(stderr,
	  "Usage: mpi_launcher [-n|-np NPROCS] ./struct --tasks INP JNP --size NI NJ NK [options]\n"
	  "    NPROCS : # of tasks launched by MPI; may or may not be implied or required by system\n\n"
	  "  Required:\n"
	  "    --tasks INP JNP: Specifies the parallel decomposition of tasks\n"
	  "      INP : # of tasks along the I (X) axis\n"
	  "      JNP : # of tasks along the J (Y) axis\n"
	  "        NOTE that INP * JNP * KNP == NPROCS is required!\n"
	  "    --size NI NJ NK : Specifies the size of the grid\n"
	  "      NI, NJ, NK : Number of grid points along the I,J,K axes respectively\n"
	  "      valid values are > 1\n\n" 
	  "  Optional:\n"
	  "    --debug : Turns on debugging print statements \n"
	  "    --maskthreshold MT : Mask theshold; valid values are floats between -1.0 and 1.0 \n"
	  "      MT : mask threshold value; Default: 0.0\n"
	  "    --noisespacefreq FNS : Spatial frequency of noise function\n"
	  "      FNS : space frequency value; Default: 10.0\n"
	  "    --noisetimefreq FNT : Temporal frequency of noise function\n"
	  "      FNT : time frequency value;  Default: 0.25\n"
	  "    --tsteps NT : Number of time steps; valid values are > 0 (Default value 10)\n"
	  "    --tstart TS : Starting time step; valid values are >= 0  (Default value 0)\n"

#ifdef HAS_ADIOS
	  "    --adios [POSIX|MPI|MPI_LUSTRE|MPI_AGGREGATE|PHDF5]: Enable ADIOS output\n"
	  "    --debugIO : Turns on debugging IO (corrently only works with ADIOS IO) \n"
#endif

	  
#ifdef HAS_HDF5
	  "    --hdf5 : Enable HDF5 output (i.e. XDMF)\n"
	  "    --hdf5_chunk y z : Chunk Size y z \n"
		  "      valid values are  NJ/JNP/y,NK/z\n"
	  "    --hdf5_compress : enable compression \n"
#endif
	  );
}
