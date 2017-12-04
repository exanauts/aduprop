#ifndef ADUPROP_ALG_HPP_
#define ADUPROP_ALG_HPP_
/*!
   \file "alg.hpp"
   \brief "Linear algebra module"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <codi.hpp>
#include "linsolve.hpp"

// First, second and third order AD datatype
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

// Passive vector definition

class pVector {
 public:
  explicit pVector(const size_t nvals);
  ~pVector();

  void set(const size_t i, const double val);
  double get(const size_t i);

  void display();

 private:
  size_t n;
  double* data;
};

// Passive matrix definition

class pMatrix {
 public:
  explicit pMatrix(const size_t nrows, const size_t ncols);
  ~pMatrix();

  void set(const size_t i, const size_t j, const double val);
  void set_col(const size_t j, const double *vals);
  double get(const size_t i, const size_t j);
  void display();

 private:
  size_t rows, cols;
  double* data;
};

#endif  // ADUPROP_ALG_HPP_
