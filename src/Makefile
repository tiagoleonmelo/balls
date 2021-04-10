CFLAGS = -g -Wall -fopenmp -O3

ALL: ballAlg ballAlg-omp ballQuery

ballAlg: ballAlg.c
	gcc $(CFLAGS) $< -lm -o $@

ballAlg-omp: ballAlg-omp.c
	kinst-ompp gcc $(CFLAGS) $< -lm -o $@

ballQuery: ballQuery.c
	gcc $< -lm -o $@

clean:
	rm -f ballAlg ballAlg-omp ballQuery