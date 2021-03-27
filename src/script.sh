#!/bin/bash

make
./"$1" 2 5 0 > output_script.txt
diff output_script.txt output_new.txt

./"$1" 2 8 0 > ex-2-8-0.tree
./ballQuery ex-2-8-0.tree 8 8 > ex-2-8-0_new.query
diff ex-2-8-0_new.query ex-2-8-0.query

set -v #echo on
./"$1" 50 1000000 0 > /dev/null
./"$1" 3 5000000 0 > /dev/null