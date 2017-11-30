#include "alg.hpp"
#include "linsolve.hpp"
#include <iostream>
#include <iomanip>

pVector::pVector(size_t nvals) {
  data = (double*)malloc(nvals*sizeof(double));
  n = nvals;
}

pVector::~pVector() {
  if (data) free(data);
}


pMatrix::pMatrix(size_t nrows, size_t ncols) {
  data = (double*)malloc(ncols*nrows*sizeof(double));
  cols = ncols;
  rows = nrows;
}

pMatrix::~pMatrix() {
  if (data) free(data);
}

void pMatrix::set(size_t i, size_t j, double val) {
  data[j*rows + i] = val;
}

double pMatrix::get(size_t i, size_t j) {
  return data[j*rows + i];
}

void pMatrix::display() {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std::cout << std::setw(5) << std::setfill(' ') << std::setprecision(2)
                  << get(i, j);
    }
    std::cout << "\n";
  }
}
