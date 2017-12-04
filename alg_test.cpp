#include <math.h>
#include <iostream>
#include "alg.hpp"


int main(int nargs, char** args) {
  
  size_t dim = 4;

  //  Vector creation
  alg::pVector b(dim);

  for (size_t i = 0; i < dim; ++i) {
    b.set(i, 1.0);
  }

  b.display();

  // Matrix creation

  alg::pMatrix A(dim, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      A.set(i, j, static_cast<double>(i + j));
    }
  }

  A.display();

  // Decremental multiplication

  alg::pVector y(dim);
 
  y.zeros();
  
  y.display();

  alg::decmatmul(A, b, y);

  y.display();

}
