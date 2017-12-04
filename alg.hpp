#ifndef ADUPROP_ALG_HPP_
#define ADUPROP_ALG_HPP_
/*!
   \file "alg.hpp"
   \brief "Linear algebra module"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <codi.hpp>

namespace alg {

// First, second and third order AD datatype
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

// Passive vector definition

template <typename T> class pVector {
 public:
  pVector();
  explicit pVector(const size_t nvals);
  ~pVector();

  void set(const size_t i, const T val);
  T get(const size_t i) const;
  size_t dim() const;
  double* get_datap() const;
  void zeros();

  void display();

 private:
  size_t n;
  double* data;
};

template <typename T> inline size_t pVector<T>::dim() const {
  return n;
}

template <> inline double* pVector<double>::get_datap() const {
  return data;
}

template <typename T> inline void pVector<T>::zeros() {
  for (size_t i = 0; i < n; ++i) {
    data[i] = 0.0;
  }
}

// Passive matrix definition

template <typename T> class pMatrix {
 public:
  pMatrix();
  explicit pMatrix(const size_t nrows, const size_t ncols);
  ~pMatrix();

  void set(const size_t i, const size_t j, const T val);
  void set_col(const size_t j, const T *vals);
  double get(const size_t i, const size_t j);
  double* get_datap() const;
  size_t nrows() const;
  size_t ncols() const;
  void display();

 private:
  size_t rows, cols;
  T* data;
};

template <typename T> inline size_t pMatrix<T>::nrows() const {
  return rows;
}

template <typename T> inline size_t pMatrix<T>::ncols() const {
  return cols;
}

template <typename T> inline double* pMatrix<T>::get_datap() const {
  return data;
}

// Function declarations

void decmatmul(const pMatrix<double> &A, const pVector<double> &x, pVector<double> &y);
void LUsolve(pMatrix<double> &A, pVector<double> &b);

} // end of namespace

#endif  // ADUPROP_ALG_HPP_
