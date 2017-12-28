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

struct Load {
  int bus;
  double P;
  double Q;
  
  Load():bus(0), P(0.0), Q(0.0){};
  void set(int bus_, double P_, double Q_);
};

void Load::set(int bus_, double P_, double Q_) {
  bus = bus_;
  P = P_;
  Q = Q_;
}

struct Generator {

  int bus;

  double x_d;
  double x_q;
  double x_dp;
  double x_qp;
  double x_ddp;
  double x_qdp;
  double xl;
  double H;
  double T_d0p;
  double T_q0p;
  double T_d0dp;
  double T_q0dp;

  double e_fd;
  double p_m;


  Generator();
  void set(int busn, double x_d_, double x_q_, double x_dp_, double x_qp_,
      double x_ddp_, double x_qdp_, double xl_, double H_, double T_d0p_,
      double T_q0p_, double T_d0dp_, double T_q0dp_);
};

Generator::Generator() {
}

void Generator::set(int busn, double x_d_, double x_q_, double x_dp_,
    double x_qp_, double x_ddp_, double x_qdp_, double xl_,
    double H_, double T_d0p_, double T_q0p_, double T_d0dp_,
    double T_q0dp_) {
  bus = busn;
  x_d = x_d_;  
  x_q = x_q_;  
  x_dp = x_dp_;  
  x_qp = x_qp_;  
  x_ddp = x_ddp_;  
  x_qdp = x_qdp_;  
  xl = xl_;  
  H = H_;  
  T_d0p = T_d0p_;  
  T_q0p = T_q0p_;  
  T_d0dp = T_d0dp_;
  T_q0dp = T_q0dp_;

  e_fd = -1.0;
  p_m = -1.0;
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
  int ngens;
  int nloads;

  double h; 
  size_t dimension;
  
  Branch* branches;
  Generator* gens;
  Load* loads;
  
  alg::pMatrix<double>* ybus;

  // Constructors
  System();
  System(int nbuses, int nbranches, int ngens, int nloads);
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
  ngens = 0;
  nloads = 0;
}


System::System(int _nbuses, int _nbranches, int _ngens, int _nloads) {
  dimension = 0; // initialize to 0
  nbuses = _nbuses;
  nbranches = _nbranches;
  ngens = _ngens;
  nloads = _nloads;

  branches = new Branch[nbranches];
  gens = new Generator[ngens];
  loads = new Load[nloads];

  // Calculate dimension
  dimension += nbuses*2; // transmission system

}

System::~System(){
  if (nbranches) delete branches;
  if (ngens) delete gens;
  if (nloads) delete loads;
  if (ybus) delete ybus;
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
