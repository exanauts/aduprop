#include "alg.hpp"
#include "linsolve.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>


namespace alg {

// VECTOR FUNCTION DEFINITIONS

pVector::pVector() {
  data = NULL;
  n = 0;
}

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

pMatrix::pMatrix() {
  data = NULL;
  cols = 0;
  rows = 0;
}

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



// LINEAR ALGEBRA FUNCTIONS

void decmatmul(const pMatrix &A, const pVector &x, pVector &y) {
  // check for dimension consistency
  assert(A.ncols() == x.dim());
  assert(A.nrows() == y.dim());
  
  // BLAS option flags
  int n = A.ncols();
  int m = A.nrows();
  char trans = 'N';
  double alpha = -1.0;
  double beta = 1.0;
  int incx = 1;
  int incy = 1;

  double* A_data = A.get_datap();
  double* x_data = x.get_datap();
  double* y_data = y.get_datap();

  FNAME(dgemv)(&trans, &m, &n, &alpha, A_data, &n, x_data,
      &incx, &beta, y_data, &incy);

}

void LUsolve(pMatrix &A, pVector &b) {
  // check for dimension consistency
  assert(A.ncols() == b.dim());
  assert(A.ncols() == A.nrows()); // square system
  
  // BLAS option flags
  int n = A.ncols();
  double* A_data = A.get_datap();
  double* b_data = b.get_datap();
  int nrhs = 1;
  int info = 0;
  int* ipiv = reinterpret_cast<int*>(calloc(n, sizeof(int)));
  int lda = n;
  int ldb = n;

  FNAME(dgesv)(&n, &nrhs, A_data, &lda, ipiv, b_data, &ldb, &info);
  assert(!info);
}

} //  End of namespace

//
