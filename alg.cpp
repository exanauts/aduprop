#include "alg.hpp"
#include <iostream>
#include <iomanip>

// VECTOR FUNCTION DEFINITIONS

pVector::pVector(const size_t nvals) {
  data = reinterpret_cast<double*>(malloc(nvals*sizeof(double)));
  n = nvals;
}

pVector::~pVector() {
  if (data) free(data);
}

void pVector::set(const size_t i, const double val) {
  data[i] = val;
}

double pVector::get(const size_t i) {
  return data[i];
}

void pVector::display() {
  for (size_t i = 0; i < n; ++i) {
    std::cout << get(i) << std::endl;
  }
}

// MATRIX DEFINITIONS

pMatrix::pMatrix(const size_t nrows, const size_t ncols) {
  data = reinterpret_cast<double*>(malloc(ncols*nrows*sizeof(double)));
  cols = ncols;
  rows = nrows;
}

pMatrix::~pMatrix() {
  if (data) free(data);
}

void pMatrix::set(const size_t i, const size_t j, const double val) {
  data[j*rows + i] = val;
}

double pMatrix::get(const size_t i, const size_t j) {
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
