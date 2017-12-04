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

template <class T> int adlinsolve(T **A, T *B, size_t n);

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
template <class T> void residual_beuler(const alg::pVector<T> &x, const alg::pVector<T> &xold,
    const double h, alg::pVector<T> &F) {

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

  const T e_qp     = x.get(0);
  const T e_dp     = x.get(1);
  const T phi_1d   = x.get(2);
  const T phi_2q   = x.get(3);
  const T w        = x.get(4);
  const T delta    = x.get(5);
  const T v_q      = x.get(6);
  const T v_d      = x.get(7);
  const T i_q      = x.get(8);
  const T i_d      = x.get(9);
  const T v1m      = x.get(10);
  const T v1a      = x.get(11);

  // Auxiliary variables.
  T psi_de, psi_qe;

  psi_de = (x_ddp - xl)/(x_dp - xl)*e_qp + 
    (x_dp - x_ddp)/(x_dp - xl)*phi_1d;

  psi_qe = -(x_ddp - xl)/(x_qp - xl)*e_dp + 
    (x_qp - x_ddp)/(x_qp - xl)*phi_2q;

  // Machine states
  F.set(0, (-e_qp + e_fd - (i_d - (-x_ddp + x_dp)*(-e_qp + i_d*(x_dp - xl) 
    + phi_1d)/pow((x_dp - xl), 2.0))*(x_d - x_dp))/T_d0p);
  F.set(1, (-e_dp + (i_q - (-x_qdp + x_qp)*( e_dp + i_q*(x_qp - xl) 
    + phi_2q)/pow((x_qp - xl), 2.0))*(x_q - x_qp))/T_q0p);
  F.set(2, ( e_qp - i_d*(x_dp - xl) - phi_1d)/T_d0dp);
  F.set(3, (-e_dp - i_q*(x_qp - xl) - phi_2q)/T_q0dp);
  F.set(4, (p_m - psi_de*i_q + psi_qe*i_d)/(2.0*H));
  F.set(5, 2.0*M_PI*60.0*w);

  // Stator currents
  F.set(6, i_d - ((x_ddp - xl)/(x_dp - xl)*e_qp + 
    (x_dp - x_ddp)/(x_dp - xl)*phi_1d - v_q)/x_ddp);
  F.set(7, i_q - (-(x_qdp - xl)/(x_qp - xl)*e_dp + 
    (x_qp - x_qdp)/(x_qp - xl)*phi_2q + v_d)/x_qdp);

  // Stator voltages
  F.set(8, v_d - v1m*sin(delta - v1a));
  F.set(9, v_q - v1m*cos(delta - v1a));

  F.set(10, v_d*i_d + v_q*i_q - ((v1m*v0m)/xline)*sin(v1a - v0a));
  F.set(11, v_q*i_d - v_d*i_q - (v1m*v1m)/xline + ((v1m*v0m)/xline)*cos(v1a - v0a));

  // Update residual for ODE equations
  F.set(0, x.get(0) - xold.get(0) - h*F.get(0));
  F.set(1, x.get(1) - xold.get(1) - h*F.get(1));
  F.set(2, x.get(2) - xold.get(2) - h*F.get(2));
  F.set(3, x.get(3) - xold.get(3) - h*F.get(3));
  F.set(4, x.get(4) - xold.get(4) - h*F.get(4));
  F.set(5, x.get(5) - xold.get(5) - h*F.get(5));

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
    const alg::pVector<T> &xold, const double h, alg::pMatrix<T> &J) {

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

  const T e_qp     = x.get(0);
  const T e_dp     = x.get(1);
  const T phi_1d   = x.get(2);
  const T phi_2q   = x.get(3);
  const T w        = x.get(4);
  const T delta    = x.get(5);
  const T v_q      = x.get(6);
  const T v_d      = x.get(7);
  const T i_q      = x.get(8);
  const T i_d      = x.get(9);
  const T v1m      = x.get(10);
  const T v1a      = x.get(11);


  // Auxiliary variables.
  T psi_de, psi_qe;


  // auxiliary variables
  psi_de = (x_ddp - xl)/(x_dp - xl)*e_qp + 
      (x_dp - x_ddp)/(x_dp - xl)*phi_1d;

  psi_qe = -(x_ddp - xl)/(x_qp - xl)*e_dp + 
      (x_qp - x_ddp)/(x_qp - xl)*phi_2q;

  // Machine states

  J.set(0, 0, 1.0 - h*(-(x_d - x_dp)*(-x_ddp + x_dp)*pow(x_dp - xl, -2.0) - 1.0)/T_d0p);
  J.set(2, 0, -h*(x_d - x_dp)*(-x_ddp + x_dp)*pow(x_dp - xl, -2.0)/T_d0p);
  J.set(9, 0, h*(x_d - x_dp)*(-(-x_ddp + x_dp)*pow(x_dp - xl, -1.0) + 1)/T_d0p );

  J.set(1, 1, 1.0 - h*(-(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0) - 1.0)/T_q0p );
  J.set(3, 1, h*(x_q - x_qp)*(-x_qdp + x_qp)*pow(x_qp - xl, -2.0)/T_q0p );
  J.set(8, 1, -h*(x_q - x_qp)*(-(-x_qdp + x_qp)*pow(x_qp - xl, -1.0) + 1.0)/T_q0p );
  
  J.set(0, 2, -h/T_d0dp );
  J.set(2, 2, 1.0 + h/T_d0dp );
  J.set(9, 2, -h*(-x_dp + xl)/T_d0dp );
  
  J.set(1, 3, h/T_q0dp );
  J.set(3, 3, 1.0 + h/T_q0dp );
  J.set(8, 3, -h*(-x_qp + xl)/T_q0dp );
  
  J.set(0, 4, 0.5*h*i_q*(x_ddp - xl)/(H*(x_dp - xl)) );
  J.set(1, 4, -0.5*h*i_d*(-x_ddp + xl)/(H*(x_qp - xl)) );
  J.set(2, 4, 0.5*h*i_q*(-x_ddp + x_dp)/(H*(x_dp - xl)) );
  J.set(3, 4, -0.5*h*i_d*(-x_ddp + x_qp)/(H*(x_qp - xl)) );
  J.set(4, 4, 1.0 );
  J.set(8, 4, -0.5*h*(-e_qp*(x_ddp - xl)/(x_dp - xl) - phi_1d*(-x_ddp + x_dp)/(x_dp - xl))/H );
  J.set(9, 4, -0.5*h*(e_dp*(-x_ddp + xl)/(x_qp - xl) + phi_2q*(-x_ddp + x_qp)/(x_qp - xl))/H) ;
   
  J.set(4, 5, -120.0*M_PI*h );
  J.set(5, 5, 1.0 );
   
  J.set(0, 6, -(x_ddp - xl)/(x_ddp*(x_dp - xl)) );
  J.set(2, 6, -(-x_ddp + x_dp)/(x_ddp*(x_dp - xl)) );
  J.set(6, 6, 1.0/x_ddp );
  J.set(9, 6, 1.0 );
   
  J.set(1, 7, -(-x_qdp + xl)/(x_qdp*(x_qp - xl)) );
  J.set(3, 7, -(-x_qdp + x_qp)/(x_qdp*(x_qp - xl)) );
  J.set(7, 7, -1/x_qdp );
  J.set(8, 7, 1.0 );
   
  J.set(5, 8, -v1m*cos(delta - v1a) );
  J.set(7, 8, 1.0 );
  J.set(10, 8, -sin(delta - v1a) );
  J.set(11, 8,  v1m*cos(delta - v1a) );
   
  J.set(5, 9,  v1m*sin(delta - v1a) );
  J.set(6, 9 , 1.0 );
  J.set(10, 9,  -cos(delta - v1a) );
  J.set(11, 9,  -v1m*sin(delta - v1a) );
  
  J.set(6, 10, i_q );
  J.set(7, 10,  i_d );
  J.set(8, 10,  v_q );
  J.set(9, 10,  v_d );
  J.set(10, 10,  v0m*sin(v0a - v1a)/xline );
  J.set(11, 10,  -v0m*v1m*cos(v0a - v1a)/xline );
  
  J.set(6, 11,  i_d );
  J.set(7, 11,  -i_q );
  J.set(8, 11,  -v_d );
  J.set(9, 11,  v_q );
  J.set(10, 11,  v0m*cos(v0a - v1a)/xline - 2.0*v1m/xline );
  J.set(11, 11,  v0m*v1m*sin(v0a - v1a)/xline );

}

/*!
   \brief "User provided time integration"
   \param x "1st order active input variable"
   \param h "Discretization"
   \pre "Initial conditions with input tangents"
   \post "New state x with 1st order tangents"
*/
template <class T> void integrate(T *x, size_t dim, double h) {
  
  double eps = 1e-9;
  int iteration = 0, ierr;
  T *xold = new T[dim];
  T *y = new T[dim];
  T **J = new T* [dim];
  
  for (size_t i = 0; i < dim; ++i)
    xold[i] = x[i];
  
  residual_beuler<T>(x, xold, h, y);
  
  J[0] = new T [dim*dim];
  
  for(size_t i = 0; i < dim; ++i)
    J[i] = J[0] + dim * i;

  for(size_t i = 0; i < dim*dim; ++i)
    J[0][i] = 0;
  
  jac_beuler<T>(x, xold, h, J);
  
  ierr = adlinsolve<T>(J, y, dim);
  
  if(ierr) {
    cout << "Linear solver error: " << ierr << endl;
    exit(1);
  }

  for(size_t i = 0; i < dim; ++i) {
    x[i] = x[i] - y[i];
  }

  delete [] xold;
  delete [] y;
  delete [] J[0];
  delete [] J;
  
} 
#endif
