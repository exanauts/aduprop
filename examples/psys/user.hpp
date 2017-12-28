/*!
  \file "user.hpp"
   \brief "This is the provided user code."
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/
#ifndef USER_HPP
#define USER_HPP
#include <iostream>
#include "alg.hpp"


// System Data structures

struct Branch {
  int fr;
  int to;
  double r;
  double x;

  Branch():fr(0), to(0), r(0), x(0){};
  Branch(int fr, int to, double r, double x):fr(fr), to(to), r(r), x(x){}
  void set(int fr_, int to_, double r_, double x_);
};

void Branch::set(int fr_, int to_, double r_, double x_) {
  fr = fr_;
  to = to_;
  r = r_;
  x = x_;
}

// Auxiliary functions

void expandComplex(int i, int j, double a, double b,
    alg::pMatrix<double>& mat) {
  mat[2*i][2*j] += a;
  mat[2*i + 1][2*j + 1] += a;
  mat[2*i][2*j + 1] += -b;
  mat[2*i + 1][2*j] += b;
}


class System {
public:
  int nbuses;
  int nbranches;
  int ngen;

  double h; 
  size_t dimension;
  
  Branch* branches;
  alg::pMatrix<double>* ybus;

  // Constructors
  System();
  System(int nbuses, int nbranches, int ngens);
  ~System();

  // Electrical system
  void build_ybus();


size_t dim() { return dimension; }

/*!
   \brief "Generic residual function written by the user. All variables that
   are on the computational path from an independent to dependent need to be of
   type T. All other variables may be passive (double)."
   \param x "New x"
   \param xold "Old x"
   \param F "Residual"
   \param h "Discretization step"
   \pre "System new state x and old state xold"
   \post "Residual F"
*/
template <class T> void residual_beuler(const alg::pVector<T> &x,
    const alg::pVector<T> &xold, alg::pVector<T> &F) {
  
  double yre, yim;

  F.zeros();

  std::cout << "Cacota" << std::endl;

  for (size_t i = 0; i < nbuses; ++i) {
    for (size_t j = 0; j < nbuses; ++j) {

      yre = (*ybus)[2*i][2*j];
      yim = (*ybus)[2*i + 1][2*j];

      if (i == j) {
        F[2*i] += x[2*i]*x[2*i]*yre;
        F[2*i + 1] += -x[2*i]*x[2*i]*yim;
      } else {
        F[2*i] += x[2*i]*x[2*j]*(yim*sin(x[2*i + 1] -
              x[2*j + 1]) + yre*cos(x[2*i + 1] - x[2*j + 1]));
        F[2*i + 1] += x[2*i]*x[2*j]*(yre*sin(x[2*i + 1] -
              x[2*j + 1]) - yim*cos(x[2*i + 1] - x[2*j + 1]));
      }
    }
  }

}

/*!
   \brief "User provided Jacobian"
   \param x "New x"
   \param xold "Old x"
   \param J "Jacobian"
   \param h "Discretization step"
   \pre "Input states x and xold"
   \post "Jacobian J"
*/
template <class T> void jac_beuler(const alg::pVector<T> &x, 
    const alg::pVector<T> &xold, alg::pMatrix<T> &J) {
  // Empty
}

};

System::System() {
  dimension = 0;
  nbranches = 0;
  nbuses = 0;
  ngen = 0;
}


System::System(int _nbuses, int _nbranches, int _ngens) {
  dimension = 0; // initialize to 0
  nbuses = _nbuses;
  nbranches = _nbranches;
  ngen = _ngens;

  branches = new Branch[nbranches];

  // Calculate dimension
  dimension += nbuses*2; // transmission system

}

System::~System(){
  if (nbranches)
    delete branches;
  if (ybus)
    delete ybus;
}


void System::build_ybus() {
  double yre, yim, mag;
  int fr, to;

  ybus = new alg::pMatrix<double>(2*nbuses, 2*nbuses);
  ybus->zeros();
  for (size_t i = 0; i < nbranches; ++i) {
    mag = pow(branches[i].r, 2.0) + pow(branches[i].x, 2.0);
    fr = branches[i].fr;
    to = branches[i].to;
    yre = branches[i].r / mag;
    yim = -branches[i].x / mag;

    expandComplex(fr, fr, yre, yim, *ybus);
    expandComplex(to, to, yre, yim, *ybus);
    expandComplex(fr, to, -yre, -yim, *ybus);
    expandComplex(to, fr, -yre, -yim, *ybus);
  }
}


#endif
