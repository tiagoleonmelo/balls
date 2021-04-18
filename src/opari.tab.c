#include "pomp_lib.h"


extern struct ompregdescr omp_rd_1;
extern struct ompregdescr omp_rd_2;
extern struct ompregdescr omp_rd_3;

int POMP_MAX_ID = 4;

struct ompregdescr* pomp_rd_table[4] = {
  0,
  &omp_rd_1,
  &omp_rd_2,
  &omp_rd_3,
};
