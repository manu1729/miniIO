mpicc -o piov1 utils.c pio-many-files-per-rank.c
mpicc -o piov2 utils.c pio-one-file-per-rank.c
mpicc -o piov3 utils.c pio-one-file-per-rank-single-write.c
