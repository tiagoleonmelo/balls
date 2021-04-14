#!/bin/bash

make
OMP_NUM_THREADS="$2" OMP_NESTED=true ./"$1" 2 5 0 > output_script.txt
diff output_script.txt output_new.txt

OMP_NUM_THREADS="$2" OMP_NESTED=true ./"$1" 2 8 0 > ex-2-8-0.tree
./ballQuery ex-2-8-0.tree 8 8 > ex-2-8-0_new.query
diff ex-2-8-0_new.query ex-2-8-0.query

set -v #echo on
OMP_NUM_THREADS="$2" OMP_NESTED=true ./"$1" 50 1000000 0 > /dev/null
#OMP_NUM_THREADS="$2" OMP_NESTED=true ./"$1" 3 5000000 0 > /dev/null
OMP_NUM_THREADS="$2" OMP_NESTED=true ./"$1" 4 20000000 0 > /dev/null