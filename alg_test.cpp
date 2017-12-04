#include <math.h>
#include <iostream>
#include "alg.hpp"


int main(int nargs, char** args) {
  
  size_t dim = 10;

  //  Vector creation
  pVector b(dim);

  for (size_t i = 0; i < dim; ++i) {
    b.set(i, 1.0);
  }

  b.display();

  // Matrix creation

  pMatrix A(dim, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      A.set(i, j, static_cast<double>(i + 10*j));
    }
  }

  A.display();
}
