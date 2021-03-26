CFLAGS = -g -Wall -fopenmp -O3

ALL: ballAlg ballQuery

ballAlg: ballAlg.c
	gcc $(CFLAGS) $< -lm -o $@

ballQuery: ballQuery.c
	gcc $< -lm -o $@

clean:
	rm -f ballAlg ballQuery