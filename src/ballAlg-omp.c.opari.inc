#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallel", "", 0, "ballAlg-omp.c", 475, 475, 479, 479
};

#define POMP_DLIST_00001 shared(omp_rd_1,omp_rd_2)

struct ompregdescr omp_rd_2 = {
  "single", "", 0, "ballAlg-omp.c", 476, 476, 479, 479
};

