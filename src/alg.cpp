#include "linsolve.hpp"
#include "alg.hpp"
#include <iostream>
#include <iomanip>
#include <cassert>
#ifdef HDF5
#include "hdf5.h"
#include "hdf5_hl.h"
#endif


namespace alg {

// VECTOR FUNCTION DEFINITIONS

template <typename T> pVector<T>::pVector() {
  data = NULL;
  n = 0;
}

template <typename T> pVector<T>::pVector(const size_t nvals) {
  data = new T[nvals];
  n = nvals;
}

template <typename T> pVector<T>::~pVector() {
  if (data != NULL) {
    delete [] data;
    data = NULL;
  } 
}

template <typename T> void pVector<T>::alloc(const size_t nvals) {
  assert(n == 0);
  assert(data == NULL);
  data = new T[nvals];
  n = nvals;
}

template <typename T> void pVector<T>::set(const size_t i, const T val) {
  data[i] = val;
}

template <typename T> T pVector<T>::get(const size_t i) const {
  return data[i];
}

#ifdef HDF5
template <> void pVector<double>::to_hdf5(const std::string filename) {
  hid_t file_id;
  hsize_t dims[2];
  dims[0] = n;
  file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset(file_id, "vector", 1, dims, H5T_NATIVE_DOUBLE, data);
  H5Fclose (file_id);
}

template <> void pVector<double>::from_hdf5(const std::string filename) {
  assert(n == 0); // vector must be empty
  hid_t file_id;
  hsize_t dims[2];
  file_id = H5Fopen (filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  H5LTget_dataset_info(file_id, "vector", dims, NULL, NULL);
  n = dims[0];
  data = (double*)malloc(dims[0]*sizeof(double));
  H5LTread_dataset_double(file_id, "vector", data);
  H5Fclose(file_id);
}
#endif

// MATRIX DEFINITIONS

template <typename T> pMatrix<T>::pMatrix() {
  data = NULL;
  cols = 0;
  rows = 0;
}

template <typename T> pMatrix<T>::pMatrix(const size_t nrows, const size_t ncols) {
  data = new T[ncols*nrows];
  cols = ncols;
  rows = nrows;
}

template <typename T> void pMatrix<T>::alloc(const size_t nrows, const size_t ncols) {
  assert(rows == 0);
  assert(cols == 0);
  assert(data == NULL);
  data = new T[ncols*nrows];
  cols = ncols;
  rows = nrows;
}

template <typename T> pMatrix<T>::~pMatrix() {
  if (data != NULL) {
    delete [] data;
    data = NULL;
  }
}

template <typename T> void pMatrix<T>::set(const size_t i, const size_t j, const T val) {
  data[j*rows + i] = val;
}


template <typename T> void pMatrix<T>::display() {
  for (size_t i = 0; i < rows; ++i) {
    for (size_t j = 0; j < cols; ++j) {
      std::cout << get(i, j) << " ";
    }
    std::cout << "\n";
  }
}


#ifdef HDF5
template <> void pMatrix<double>::to_hdf5(const std::string filename) {
  hid_t file_id;
  hsize_t dims[2];
  dims[0] = rows;
  dims[1] = cols;
  file_id = H5Fcreate (filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  H5LTmake_dataset(file_id, "matrix", 2, dims, H5T_NATIVE_DOUBLE, data);
  H5Fclose (file_id);
}
#endif

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
  bool first = true;
  int* ipiv = NULL;
  if(first) {  
    ipiv = reinterpret_cast<int*>(calloc(n, sizeof(int)));
    first = false;
  }
  int lda = n;
  int ldb = n;

  FNAME(dgesv)(&n, &nrhs, A_data, &lda, ipiv, b_data, &ldb, &info);
  free(ipiv);
  assert(!info);
}

template class pVector<double>;
template class pVector<t1s>;
template class pVector<t2s>;
template class pVector<t3s>;
template class pMatrix<double>;
template class pMatrix<t1s>;
template class pMatrix<t2s>;
template class pMatrix<t3s>;

} //  End of namespace

//
