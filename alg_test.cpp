#include <math.h>
#include <iostream>
#include "alg.hpp"

#include <cassert>

int main(int nargs, char** args) {
  
  size_t dim = 4;

  double EPS = 1e-10;

  //  Vector creation
  alg::pVector<double> b(dim);

  for (size_t i = 0; i < dim; ++i) {
    b.set(i, 1.0);
  }

  std::cout << b << std::endl;

  // Matrix creation

  alg::pMatrix<double> A(dim, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      A.set(i, j, static_cast<double>(pow(i + 1, j)));
    }
  }

  std::cout << A << std::endl;

  // Decremental multiplication

  alg::pVector<double> y(dim);
  y.zeros();
  alg::decmatmul(A, b, y);
  std::cout << y << std::endl;

  // LU decomp

  alg::LUsolve(A, b);
  std::cout << b << std::endl;


  // Performance test
  

  uint64_t cycleSTART, cycleEND;

  size_t large_dim = 10000;

  alg::pMatrix<double> B(large_dim, large_dim);
  double foo = 0;

  cycleSTART = __rdtsc();
  // Vanilla accessor
  for (size_t i = 0; i < large_dim; ++i) {
    for (size_t j = 0; j < large_dim; ++j) {
      foo += B.get(i, j);
    }
  }
  cycleEND = __rdtsc();

  std::cout << "Number of cycles spent in vanilla: " << cycleEND - cycleSTART << std::endl;

  foo = 0;

  // Operator overload
  
  cycleSTART = __rdtsc();
  for (size_t i = 0; i < large_dim; ++i) {
    for (size_t j = 0; j < large_dim; ++j) {
      foo += B[i][j];
    }
  }
  cycleEND = __rdtsc();
  
  std::cout << "Number of cycles spent in overload: " << cycleEND - cycleSTART << std::endl;


  // IO stuff
 
#ifdef HDF5
  alg::pVector<double> b2, err;
  b.to_hdf5("b.hdf5");
  b2.from_hdf5("b.hdf5");
  err = b - b2;
  assert(abs(err.norm()) < EPS);
#endif

  


}
