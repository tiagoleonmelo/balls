#!/bin/bash

for i in 1 2 4 8 16 32 64
do
   ./tablegen.sh $i ballAlg-mpi
   ./tablegen.sh $i ballAlg-mpi-openmp
   ./tablegen.sh $i ballAlg-mpi-openmpFor
   #./tablegen.sh $i ballAlg-mpi-openmpTasks
done
