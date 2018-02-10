#ifndef ADUPROP_PARALLEL_HPP_
#define ADUPROP_PARALLEL_HPP_

#include <mpi.h>

int paduprop_init();
int paduprop_destroy();
size_t paduprop_getstart(size_t dim);
size_t paduprop_getend(size_t dim);

int paduprop_getrank();
int paduprop_getcommsize();

#endif  // ADUPROP_PARALLEL_HPP_
