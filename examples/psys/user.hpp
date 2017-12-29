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


#define GEN_SIZE 10

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

  double step_size; // step size for beuler
  size_t dimension;
  size_t pnet; // pointer to network equations
  
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
  
  size_t genp;
  double yre, yim;
  T vfr, vto, afr, ato;
  T vbus, abus;
  
  // Generator temporary state vars
  T e_qp, e_dp, phi_1d, phi_2q, w, delta;
  T v_q, v_d, i_q, i_d;
  T psi_de, psi_qe;
  double x_d, x_q, x_dp, x_qp, x_ddp, x_qdp;
  double xl, H, T_d0p, T_q0p, T_d0dp, T_q0dp;
  double e_fd, p_m;

  F.zeros();

  for (size_t i = 0; i < ngens; ++i) {
    genp = i*GEN_SIZE;

    // retrieve state var values
    e_qp     = x[genp];
    e_dp     = x[genp + 1];
    phi_1d   = x[genp + 2];
    phi_2q   = x[genp + 3];
    w        = x[genp + 4];
    delta    = x[genp + 5];
    v_q      = x[genp + 6];
    v_d      = x[genp + 7];
    i_q      = x[genp + 8];
    i_d      = x[genp + 9];

    // retrieve generator parameters
    x_d = gens[i].x_d;
    x_q = gens[i].x_q;
    x_dp = gens[i].x_dp;
    x_qp = gens[i].x_qp;
    x_ddp = gens[i].x_ddp;
    x_qdp = gens[i].x_qdp;
    xl = gens[i].xl;
    H = gens[i].H;
    T_d0p = gens[i].T_d0p;
    T_q0p = gens[i].T_q0p;
    T_d0dp = gens[i].T_d0dp;
    T_q0dp = gens[i].T_q0dp;
    e_fd = gens[i].e_fd;
    p_m = gens[i].p_m;

    // Voltage
    vbus = x[pnet + 2*gens[i].bus];
    abus = x[pnet + 2*gens[i].bus + 1];

    psi_de = (x_ddp - xl)/(x_dp - xl)*e_qp + 
      (x_dp - x_ddp)/(x_dp - xl)*phi_1d;

    psi_qe = -(x_ddp - xl)/(x_qp - xl)*e_dp + 
      (x_qp - x_ddp)/(x_qp - xl)*phi_2q;

    // Machine states
    F[genp] = (-e_qp + e_fd - (i_d - (-x_ddp + x_dp)*(-e_qp + i_d*(x_dp - xl) 
      + phi_1d)/pow((x_dp - xl), 2.0))*(x_d - x_dp))/T_d0p;
    F[genp + 1] = (-e_dp + (i_q - (-x_qdp + x_qp)*( e_dp + i_q*(x_qp - xl) 
      + phi_2q)/pow((x_qp - xl), 2.0))*(x_q - x_qp))/T_q0p;
    F[genp + 2] = ( e_qp - i_d*(x_dp - xl) - phi_1d)/T_d0dp;
    F[genp + 3] = (-e_dp - i_q*(x_qp - xl) - phi_2q)/T_q0dp;
    F[genp + 4] = (p_m - psi_de*i_q + psi_qe*i_d)/(2.0*H);
    F[genp + 5] = 2.0*M_PI*60.0*w;

    // Stator currents
    F[genp + 6] = i_d - ((x_ddp - xl)/(x_dp - xl)*e_qp + 
      (x_dp - x_ddp)/(x_dp - xl)*phi_1d - v_q)/x_ddp;
    F[genp + 7] = i_q - (-(x_qdp - xl)/(x_qp - xl)*e_dp + 
      (x_qp - x_qdp)/(x_qp - xl)*phi_2q + v_d)/x_qdp;

    // Stator voltages
    F[genp + 8] = v_d - vbus*sin(delta - abus);
    F[genp + 9] = v_q - vbus*cos(delta - abus);

    F[pnet + 2*gens[i].bus] += v_d*i_d + v_q*i_q;
    F[pnet + 2*gens[i].bus + 1] += v_q*i_d - v_d*i_q;
  }



  // Power flow equations
  for (size_t i = 0; i < nbuses; ++i) {
    for (size_t j = 0; j < nbuses; ++j) {

      yre = (*ybus)[2*i][2*j];
      yim = (*ybus)[2*i + 1][2*j];

      vfr = x[pnet + 2*i];
      vto = x[pnet + 2*j];
      afr = x[pnet + 2*i + 1];
      ato = x[pnet + 2*j + 1];

      if (i == j) {
        F[pnet + 2*i] -= vfr*vto*yre;
        F[pnet + 2*i + 1] -= -vfr*vto*yim;
      } else {
        F[pnet + 2*i] -= vfr*vto*(yim*sin(afr - ato) 
            + yre*cos(afr - ato));
        F[pnet + 2*i + 1] -= vfr*vto*(yre*sin(afr - ato)
            - yim*cos(afr - ato));
      }
    }
  }

  // Load injections (PQ)
  for (size_t i = 0; i < nloads; ++i) {
    F[pnet + 2*loads[i].bus] -= loads[i].P;
    F[pnet + 2*loads[i].bus + 1] -= loads[i].Q;
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
  pnet = ngens*GEN_SIZE; // assume all gens are GENROU
  dimension = pnet + nbuses*2; // add transmission system

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
