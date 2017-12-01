#include <math.h>
#include <iostream>
#include "alg.hpp"


int main(int nargs, char** args) {
  pMatrix A(5, 5);

  for (size_t i = 0; i < 5; ++i) {
    for (size_t j = 0; j < 5; ++j) {
      A.set(i, j, static_cast<double>(i + 10*j));
    }
  }

  A.display();
  double caca = A.get(2, 2);
}
