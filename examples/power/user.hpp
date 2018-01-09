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

typedef struct System {
public:
  // System parameters.
  // Generator
  double const x_d = 1.575;
  double const x_q = 1.512;
  double const x_dp = 0.29;
  double const x_qp = 0.39;
  double const x_ddp = 0.1733;
  double const x_qdp = 0.1733;
  double const xl = 0.0787;
  double const H = 3.38;
  double const T_d0p = 6.09;
  double const T_q0p = 1.0;
  double const T_d0dp = 0.05;
  double const T_q0dp = 0.15;
  // oinfinite bus
  double const v0m = 1.0162384;
  double const v0a = -0.05807256;
  double const xline = 0.0576;
  
  const double e_fd = 2.36980307364616349375;
  const double p_m = 1.06496000000000012875;
  
  const double h = 0.004; 
  const size_t dimension = 12;

void ic(pVector<double> &x) {
  // Initial values for state array.

  x[0] = 1.06512037300928485983;
  x[1] = 0.51819992367912581788;
  x[2] = 0.85058383242985102779;
  x[3] = -0.66197500054304025952;
  x[4] = -0.05000000000000000278;
  x[5] = 0.73618306350367335167;
  x[6] = 0.77067836274882195458;
  x[7] = 0.69832289180288620312;
  x[8] = 0.46185376441989828278;
  x[9] = 1.01531727676021699125;
  x[10] = 1.0400000000000000000;
  x[11] = 0.00000000000000000000;
}
  
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
template <class T> void residual_beuler(const pVector<T> &x, const pVector<T> &xold,
    pVector<T> &F) {

  // (TEMP): Just put all the parameters here for now.

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
    const pVector<T> &xold, pMatrix<T> &J) {

  const T e_qp     = x[0];
  const T e_dp     = x[1];
  const T phi_1d   = x[2];
  const T phi_2q   = x[3];
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
  J[0][2] = -h*(x_d - x_dp)*(-x_ddp + x_dp)*pow(x_dp - xl, -2.0)/T_d0p;
  J[0][9] = h*(x_d - x_dp)*(-(-x_ddp + x_dp)*pow(x_dp - xl, -1.0) + 1)/T_d0p;

  J[1][1] = 1.0 - h*(-(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0) - 1.0)/T_q0p;
  J[1][3] = h*(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0)/T_q0p;
  J[1][8] = -h*(x_q - x_qp)*(-(-x_qdp + x_qp)*pow(x_qp - xl, -1.0) + 1.0)/T_q0p;
  
  J[2][0] = -h/T_d0dp;
  J[2][2] = 1.0 + h/T_d0dp;
  J[2][9] = -h*(-x_dp + xl)/T_d0dp;
  
  J[3][1] = h/T_q0dp;
  J[3][3] = 1.0 + h/T_q0dp;
  J[3][8] = -h*(-x_qp + xl)/T_q0dp;
  
  J[4][0] = 0.5*h*i_q*(x_ddp - xl)/(H*(x_dp - xl));
  J[4][1] = -0.5*h*i_d*(-x_ddp + xl)/(H*(x_qp - xl));
  J[4][2] = 0.5*h*i_q*(-x_ddp + x_dp)/(H*(x_dp - xl));
  J[4][3] = -0.5*h*i_d*(-x_ddp + x_qp)/(H*(x_qp - xl));
  J[4][4] = 1.0;
  J[4][8] = -0.5*h*(-e_qp*(x_ddp - xl)/(x_dp - xl) - phi_1d*(-x_ddp + x_dp)/(x_dp - xl))/H;
  J[4][9] = -0.5*h*(e_dp*(-x_ddp + xl)/(x_qp - xl) + phi_2q*(-x_ddp + x_qp)/(x_qp - xl))/H;
   
  J[5][4] = -120.0*M_PI*h;
  J[5][5] = 1.0;
   
  J[6][0] = -(x_ddp - xl)/(x_ddp*(x_dp - xl));
  J[6][2] = -(-x_ddp + x_dp)/(x_ddp*(x_dp - xl));
  J[6][6] = 1.0/x_ddp;
  J[6][9] = 1.0;
   
  J[7][1] = -(-x_qdp + xl)/(x_qdp*(x_qp - xl));
  J[7][3] = -(-x_qdp + x_qp)/(x_qdp*(x_qp - xl));
  J[7][7] = -1/x_qdp;
  J[7][8] = 1.0;
   
  J[8][5] = -v1m*cos(delta - v1a);
  J[8][7] = 1.0;
  J[8][10] = -sin(delta - v1a);
  J[8][11] = v1m*cos(delta - v1a);
   
  J[9][5] = v1m*sin(delta - v1a);
  J[9][6] = 1.0;
  J[9][10] = -cos(delta - v1a);
  J[9][11] = -v1m*sin(delta - v1a);
  
  J[10][6] = i_q;
  J[10][7] = i_d;
  J[10][8] = v_q;
  J[10][9] = v_d;
  J[10][10] = v0m*sin(v0a - v1a)/xline;
  J[10][11] = -v0m*v1m*cos(v0a - v1a)/xline;
  
  J[11][6] = i_d;
  J[11][7] = -i_q;
  J[11][8] = -v_d;
  J[11][9] = v_q;
  J[11][10] = v0m*cos(v0a - v1a)/xline - 2.0*v1m/xline;
  J[11][11] = v0m*v1m*sin(v0a - v1a)/xline;

}

} System;
#endif
