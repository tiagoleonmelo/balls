ballAlg: ballAlg.o
	gcc ballAlg.o -fopenmp -O3 -lm -o ballAlg

ballAlg.o: ballAlg.c
	gcc -c ballAlg.c -o ballAlg.o

clean:
	rm -f ballAlg.o ballAlg