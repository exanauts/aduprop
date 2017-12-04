#include "linsolve.hpp"
#include "alg.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>


namespace alg {

// VECTOR FUNCTION DEFINITIONS

template <typename T> pVector<T>::pVector() {
  data = NULL;
  n = 0;
}

template <> pVector<double>::pVector(const size_t nvals) {
  data = reinterpret_cast<double*>(malloc(nvals*sizeof(double)));
  n = nvals;
}

template <typename T> pVector<T>::pVector(const size_t nvals) {
  data = new T[nvals];
  n = nvals;
}

template <> pVector<double>::~pVector() {
  if (data) free(data);
}

template <typename T> pVector<T>::~pVector() {
  if (data) delete [] data;
}

template <typename T> void pVector<T>::set(const size_t i, const T val) {
  data[i] = val;
}

template <typename T> T pVector<T>::get(const size_t i) const {
  return data[i];
}

template <typename T> void pVector<T>::display() {
  for (size_t i = 0; i < n; ++i) {
    std::cout << get(i) << std::endl;
  }
}

// MATRIX DEFINITIONS

template <typename T> pMatrix<T>::pMatrix() {
  data = NULL;
  cols = 0;
  rows = 0;
}

template <> pMatrix<double>::pMatrix(const size_t nrows, const size_t ncols) {
  data = reinterpret_cast<double*>(malloc(ncols*nrows*sizeof(double)));
  cols = ncols;
  rows = nrows;
}

template <> pMatrix<double>::~pMatrix() {
  if (data) free(data);
}

template <typename T> void pMatrix<T>::set(const size_t i, const size_t j, const T val) {
  data[j*rows + i] = val;
}

template <typename T> T pMatrix<T>::get(const size_t i, const size_t j) {
  return data[j*rows + i];
}

template <typename T> void pMatrix<T>::display() {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std::cout << get(i, j) << " ";
    }
    std::cout << "\n";
  }
}



// LINEAR ALGEBRA FUNCTIONS

void decmatmul(const pMatrix<double> &A, const pVector<double> &x, 
    pVector<double> &y) {
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

void LUsolve(pMatrix<double> &A, pVector<double> &b) {
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

template class pVector<double>;
template class pVector<t1s>;
template class pMatrix<double>;
template class pMatrix<t1s>;

} //  End of namespace

//
