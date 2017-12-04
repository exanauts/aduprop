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

class pVector {
 public:
  pVector();
  explicit pVector(const size_t nvals);
  ~pVector();

  void set(const size_t i, const double val);
  double get(const size_t i);
  size_t dim() const;
  double* get_datap() const;
  void zeros();

  void display();

 private:
  size_t n;
  double* data;
};

inline size_t pVector::dim() const {
  return n;
}

inline double* pVector::get_datap() const {
  return data;
}

inline void pVector::zeros() {
  for (size_t i = 0; i < n; ++i) {
    data[i] = 0.0;
  }
}

// Passive matrix definition

class pMatrix {
 public:
  pMatrix();
  explicit pMatrix(const size_t nrows, const size_t ncols);
  ~pMatrix();

  void set(const size_t i, const size_t j, const double val);
  void set_col(const size_t j, const double *vals);
  double get(const size_t i, const size_t j);
  double* get_datap() const;
  size_t nrows() const;
  size_t ncols() const;
  void display();

 private:
  size_t rows, cols;
  double* data;
};

inline size_t pMatrix::nrows() const {
  return rows;
}

inline size_t pMatrix::ncols() const {
  return cols;
}

inline double* pMatrix::get_datap() const {
  return data;
}

// Function declarations

void decmatmul(const pMatrix &A, const pVector &x, pVector &y);

} // end of namespace

#endif  // ADUPROP_ALG_HPP_
