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

}
