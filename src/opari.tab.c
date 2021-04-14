#include "pomp_lib.h"


extern struct ompregdescr omp_rd_27;
extern struct ompregdescr omp_rd_28;
extern struct ompregdescr omp_rd_29;
extern struct ompregdescr omp_rd_30;

int POMP_MAX_ID = 31;

struct ompregdescr* pomp_rd_table[31] = {
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  0,
  &omp_rd_27,
  &omp_rd_28,
  &omp_rd_29,
  &omp_rd_30,
};
