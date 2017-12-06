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

using namespace std;

using namespace alg;

template <class T> void adlinsolve(pMatrix<T> &A, pVector<T> &B);

typedef struct System {
  // Generator
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
  // oinfinite bus
  double v0m;
  double v0a;
  double xline;
};

System sys;
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
template <class T> void residual_beuler(const pVector<T> &x, const pVector<T> &xold,
    const double h, pVector<T> &F) {

  // (TEMP): Just put all the parameters here for now.

  const double e_fd = 2.36980307364616349375;
  const double p_m = 1.06496000000000012875;

  const double x_d = sys.x_d; 
  const double x_q = sys.x_q;
  const double x_dp = sys.x_dp;
  const double x_qp = sys.x_qp;
  const double x_ddp = sys.x_ddp; 
  const double x_qdp = sys.x_qdp;
  const double xl = sys.xl;
  const double H = sys.H; 
  const double T_d0p = sys.T_d0p;
  const double T_q0p = sys.T_q0p;
  const double T_d0dp = sys.T_d0dp;
  const double T_q0dp = sys.T_q0dp;

  // infinite bus
  const double v0m = sys.v0m;
  const double v0a = sys.v0a;
  const double xline = sys.xline;

  const T e_qp     = x[0];
  const T e_dp     = x[1];
  const T phi_1d   = x[2];
  const T phi_2q   = x[3];
  const T w        = x[4];
  const T delta    = x[5];
  const T v_q      = x[6];
  const T v_d      = x[7];
  const T i_q      = x[8];
  const T i_d      = x[9];
  const T v1m      = x[10];
  const T v1a      = x[11];

  // Auxiliary variables.
  T psi_de, psi_qe;

  psi_de = (x_ddp - xl)/(x_dp - xl)*e_qp + 
    (x_dp - x_ddp)/(x_dp - xl)*phi_1d;

  psi_qe = -(x_ddp - xl)/(x_qp - xl)*e_dp + 
    (x_qp - x_ddp)/(x_qp - xl)*phi_2q;

  // Machine states
  F[0] = (-e_qp + e_fd - (i_d - (-x_ddp + x_dp)*(-e_qp + i_d*(x_dp - xl) 
    + phi_1d)/pow((x_dp - xl), 2.0))*(x_d - x_dp))/T_d0p;
  F[1] = (-e_dp + (i_q - (-x_qdp + x_qp)*( e_dp + i_q*(x_qp - xl) 
    + phi_2q)/pow((x_qp - xl), 2.0))*(x_q - x_qp))/T_q0p;
  F[2] = ( e_qp - i_d*(x_dp - xl) - phi_1d)/T_d0dp;
  F[3] = (-e_dp - i_q*(x_qp - xl) - phi_2q)/T_q0dp;
  F[4] = (p_m - psi_de*i_q + psi_qe*i_d)/(2.0*H);
  F[5] = 2.0*M_PI*60.0*w;

  // Stator currents
  F[6] = i_d - ((x_ddp - xl)/(x_dp - xl)*e_qp + 
    (x_dp - x_ddp)/(x_dp - xl)*phi_1d - v_q)/x_ddp;
  F[7] = i_q - (-(x_qdp - xl)/(x_qp - xl)*e_dp + 
    (x_qp - x_qdp)/(x_qp - xl)*phi_2q + v_d)/x_qdp;

  // Stator voltages
  F[8] = v_d - v1m*sin(delta - v1a);
  F[9] = v_q - v1m*cos(delta - v1a);

  F[10] = v_d*i_d + v_q*i_q - ((v1m*v0m)/xline)*sin(v1a - v0a);
  F[11] = v_q*i_d - v_d*i_q - (v1m*v1m)/xline + ((v1m*v0m)/xline)*cos(v1a - v0a);

  // Update residual for ODE equations
  F[0] = x[0] - xold[0] - h*F[0];
  F[1] = x[1] - xold[1] - h*F[1];
  F[2] = x[2] - xold[2] - h*F[2];
  F[3] = x[3] - xold[3] - h*F[3];
  F[4] = x[4] - xold[4] - h*F[4];
  F[5] = x[5] - xold[5] - h*F[5];

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
template <class T> void jac_beuler(const pVector<T> &x, 
    const pVector<T> &xold, const double h, pMatrix<T> &J) {

  size_t ndim = 12;

  const double e_fd = 2.36980307364616349375;
  const double p_m = 1.06496000000000012875;

  const double x_d = sys.x_d; 
  const double x_q = sys.x_q;
  const double x_dp = sys.x_dp;
  const double x_qp = sys.x_qp;
  const double x_ddp = sys.x_ddp; 
  const double x_qdp = sys.x_qdp;
  const double xl = sys.xl;
  const double H = sys.H; 
  const double T_d0p = sys.T_d0p;
  const double T_q0p = sys.T_q0p;
  const double T_d0dp = sys.T_d0dp;
  const double T_q0dp = sys.T_q0dp;

  // infinite bus
  const double v0m = sys.v0m;
  const double v0a = sys.v0a;
  const double xline = sys.xline;

  const T e_qp     = x[0];
  const T e_dp     = x[1];
  const T phi_1d   = x[2];
  const T phi_2q   = x[3];
  const T w        = x[4];
  const T delta    = x[5];
  const T v_q      = x[6];
  const T v_d      = x[7];
  const T i_q      = x[8];
  const T i_d      = x[9];
  const T v1m      = x[10];
  const T v1a      = x[11];


  // Auxiliary variables.
  T psi_de, psi_qe;


  // auxiliary variables
  psi_de = (x_ddp - xl)/(x_dp - xl)*e_qp + 
      (x_dp - x_ddp)/(x_dp - xl)*phi_1d;

  psi_qe = -(x_ddp - xl)/(x_qp - xl)*e_dp + 
      (x_qp - x_ddp)/(x_qp - xl)*phi_2q;

  // Machine states

  J.zeros();
  J[0][0] = 1.0 - h*(-(x_d - x_dp)*(-x_ddp + x_dp)*pow(x_dp - xl, -2.0) - 1.0)/T_d0p;
  J[2][0] = -h*(x_d - x_dp)*(-x_ddp + x_dp)*pow(x_dp - xl, -2.0)/T_d0p;
  J[9][0] = h*(x_d - x_dp)*(-(-x_ddp + x_dp)*pow(x_dp - xl, -1.0) + 1)/T_d0p;

  J[1][1] = 1.0 - h*(-(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0) - 1.0)/T_q0p;
  J[3][1] = h*(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0)/T_q0p;
  J[8][1] = -h*(x_q - x_qp)*(-(-x_qdp + x_qp)*pow(x_qp - xl, -1.0) + 1.0)/T_q0p;
  
  J[0][2] = -h/T_d0dp;
  J[2][2] = 1.0 + h/T_d0dp;
  J[9][2] = -h*(-x_dp + xl)/T_d0dp;
  
  J[1][3] = h/T_q0dp;
  J[3][3] = 1.0 + h/T_q0dp;
  J[8][3] = -h*(-x_qp + xl)/T_q0dp;
  
  J[0][4] = 0.5*h*i_q*(x_ddp - xl)/(H*(x_dp - xl));
  J[1][4] = -0.5*h*i_d*(-x_ddp + xl)/(H*(x_qp - xl));
  J[2][4] = 0.5*h*i_q*(-x_ddp + x_dp)/(H*(x_dp - xl));
  J[3][4] = -0.5*h*i_d*(-x_ddp + x_qp)/(H*(x_qp - xl));
  J[4][4] = 1.0;
  J[8][4] = -0.5*h*(-e_qp*(x_ddp - xl)/(x_dp - xl) - phi_1d*(-x_ddp + x_dp)/(x_dp - xl))/H;
  J[9][4] = -0.5*h*(e_dp*(-x_ddp + xl)/(x_qp - xl) + phi_2q*(-x_ddp + x_qp)/(x_qp - xl))/H;
   
  J[4][5] = -120.0*M_PI*h;
  J[5][5] = 1.0;
   
  J[0][6] = -(x_ddp - xl)/(x_ddp*(x_dp - xl));
  J[2][6] = -(-x_ddp + x_dp)/(x_ddp*(x_dp - xl));
  J[6][6] = 1.0/x_ddp;
  J[9][6] = 1.0;
   
  J[1][7] = -(-x_qdp + xl)/(x_qdp*(x_qp - xl));
  J[3][7] = -(-x_qdp + x_qp)/(x_qdp*(x_qp - xl));
  J[7][7] = -1/x_qdp;
  J[8][7] = 1.0;
   
  J[5][8] = -v1m*cos(delta - v1a);
  J[7][8] = 1.0;
  J[10][8] = -sin(delta - v1a);
  J[11][8] = v1m*cos(delta - v1a);
   
  J[5][9] = v1m*sin(delta - v1a);
  J[6][9] = 1.0;
  J[10][9] = -cos(delta - v1a);
  J[11][9] = -v1m*sin(delta - v1a);
  
  J[6][10] = i_q;
  J[7][10] = i_d;
  J[8][10] = v_q;
  J[9][10] = v_d;
  J[10][10] = v0m*sin(v0a - v1a)/xline;
  J[11][10] = -v0m*v1m*cos(v0a - v1a)/xline;
  
  J[6][11] = i_d;
  J[7][11] = -i_q;
  J[8][11] = -v_d;
  J[9][11] = v_q;
  J[10][11] = v0m*cos(v0a - v1a)/xline - 2.0*v1m/xline;
  J[11][11] = v0m*v1m*sin(v0a - v1a)/xline;

}

/*!
   \brief "User provided time integration"
   \param x "1st order active input variable"
   \param h "Discretization"
   \pre "Initial conditions with input tangents"
   \post "New state x with 1st order tangents"
*/
template <class T> void integrate(pVector<T> &x, size_t dim, double h) {
  
  double eps = 1e-9;
  int iteration = 0;
  pVector<T> xold(dim);
  pVector<T> y(dim);
  pMatrix<T> J(dim,dim);
  
  xold = x;
  residual_beuler<T>(x, xold, h, y);
  J.zeros();
  
  do {
    iteration = iteration + 1;
    jac_beuler<T>(x, xold, h, J);
    residual_beuler<T>(x, xold, h, y);
    adlinsolve<T>(J, y);
    x = x - y;  
    residual_beuler<T>(x, xold, h, y);
  } while (y.norm() > eps);
  
} 
#endif
