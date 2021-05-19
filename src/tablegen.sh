#!/bin/bash

filename=time"$1"_"$2".txt

set -v #echo on
srun -n "$1" ./"$2" 2 5 0 > /dev/null 2> $filename

srun -n "$1" ./"$2" 20 1000000 0 > /dev/null 2>> $filename

srun -n "$1" ./"$2" 3 5000000 0 > /dev/null 2>> $filename

srun -n "$1" ./"$2" 4 10000000 0 > /dev/null 2>> $filename

srun -n "$1" ./"$2" 3 20000000 0 > /dev/null 2>> $filename

srun -n "$1" ./"$2" 4 20000000 0 > /dev/null 2>> $filename

srun -n "$1" ./"$2" 5 20000000 0 > /dev/null 2>> $filename

#srun -n "$1" ./"$2" 6 20000000 0 > /dev/null 2>> $filename

#srun -n "$1" ./"$2" 7 20000000 0 > /dev/null 2>> $filename
set +v #echo off
