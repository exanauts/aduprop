#include <math.h>
#include <iostream>
#include "alg.hpp"

#include <cassert>

int main(int nargs, char** args) {
  
  size_t dim = 4;

  //  Vector creation
  alg::pVector<double> b(dim);

  for (size_t i = 0; i < dim; ++i) {
    b.set(i, 1.0);
  }

  b.display();

  // Matrix creation

  alg::pMatrix<double> A(dim, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      A.set(i, j, static_cast<double>(pow(i + 1, j)));
    }
  }

  A.display();

  // Decremental multiplication

  alg::pVector<double> y(dim);
  y.zeros();
  alg::decmatmul(A, b, y);
  y.display();

  // LU decomp

  alg::LUsolve(A, b);
  b.display();


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

}
