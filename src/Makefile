CFLAGS = -g -Wall -fopenmp -O3

ALL: ballAlg ballAlg-omp ballAlg-mpi

ballAlg: ballAlg.c
	gcc $(CFLAGS) $< -lm -o $@

ballAlg-omp: ballAlg-omp.c
	gcc $(CFLAGS) $< -lm -o $@

ballAlg-mpi: ballAlg-mpi.c
	mpicc $(CFLAGS) $< -lm -o $@

clean:
	rm -f ballAlg ballAlg-omp ballAlg-mpi