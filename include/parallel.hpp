#ifndef ADUPROP_PARALLEL_HPP_
#define ADUPROP_PARALLEL_HPP_

#include <mpi.h>
#include <fstream>
#include "alg.hpp"
using namespace alg;

int paduprop_init();
int paduprop_destroy();
size_t paduprop_getstart(size_t dim);
size_t paduprop_getend(size_t dim);

int paduprop_getrank();
int paduprop_getcommsize();
int paduprop_sum(pMatrix<double>&);
int paduprop_sum(pVector<double> &vec);
int paduprop_gather(pTensor3<double> &ten);
int paduprop_gather(pMatrix<double> &mat);
int paduprop_gather(pVector<double> &vec);

#endif  // ADUPROP_PARALLEL_HPP_
