#include "pomp_lib.h"

struct ompregdescr omp_rd_1 = {
  "parallelfor", "", 0, "ballAlg-omp.c", 100, 100, 104, 104
};

#define POMP_DLIST_00001 shared(omp_rd_1)

struct ompregdescr omp_rd_2 = {
  "parallelfor", "", 0, "ballAlg-omp.c", 128, 128, 146, 146
};

#define POMP_DLIST_00002 shared(omp_rd_2)

struct ompregdescr omp_rd_3 = {
  "parallelfor", "", 0, "ballAlg-omp.c", 626, 626, 631, 631
};

#define POMP_DLIST_00003 shared(omp_rd_3)

