#include <cmath>
#include <codi.hpp>
#include <iostream>
#include "linsolve.hpp"


typedef codi::RealForward active;

struct System {
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

template <class T> void residual_beuler(const T* const x, const T* const xold,
    System* const sys, const double h, T* const F) {

  // (TEMP): Just put all the parameters here for now.

  const double e_fd = 2.36980307368;
  const double p_m = 1.06496;

  const double x_d = sys->x_d; 
  const double x_q = sys->x_q;
  const double x_dp = sys->x_dp;
  const double x_qp = sys->x_qp;
  const double x_ddp = sys->x_ddp; 
  const double x_qdp = sys->x_qdp;
  const double xl = sys->xl;
  const double H = sys->H; 
  const double T_d0p = sys->T_d0p;
  const double T_q0p = sys->T_q0p;
  const double T_d0dp = sys->T_d0dp;
  const double T_q0dp = sys->T_q0dp;

  // infinite bus
  const double v0m = sys->v0m;
  const double v0a = sys->v0a;
  const double xline = sys->xline;

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
    + phi_1d)/std::pow((x_dp - xl), 2.0))*(x_d - x_dp))/T_d0p;
  F[1] = (-e_dp + (i_q - (-x_qdp + x_qp)*( e_dp + i_q*(x_qp - xl) 
    + phi_2q)/std::pow((x_qp - xl), 2.0))*(x_q - x_qp))/T_q0p;
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


template <class T> void jac_beuler(const T* const x, const T* const xold,
    System* const sys, const double h, T* const J) {
    

  size_t ndim = 12;

  const double e_fd = 2.36980307368;
  const double p_m = 1.06496;

  const double x_d = sys->x_d; 
  const double x_q = sys->x_q;
  const double x_dp = sys->x_dp;
  const double x_qp = sys->x_qp;
  const double x_ddp = sys->x_ddp; 
  const double x_qdp = sys->x_qdp;
  const double xl = sys->xl;
  const double H = sys->H; 
  const double T_d0p = sys->T_d0p;
  const double T_q0p = sys->T_q0p;
  const double T_d0dp = sys->T_d0dp;
  const double T_q0dp = sys->T_q0dp;

  // infinite bus
  const double v0m = sys->v0m;
  const double v0a = sys->v0a;
  const double xline = sys->xline;

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

  J[0*ndim + 0] = 1 - h*(-(x_d - x_dp)*(-x_ddp + x_dp)*std::pow(x_dp - xl, -2.0) - 1)/T_d0p;
  J[0*ndim + 2] = -h*(x_d - x_dp)*(-x_ddp + x_dp)*std::pow(x_dp - xl, -2.0)/T_d0p;
  J[0*ndim + 9] = h*(x_d - x_dp)*(-(-x_ddp + x_dp)*std::pow(x_dp - xl, -1.0) + 1)/T_d0p;

  J[1*ndim + 1] = 1 - h*(-(x_q - x_qp)*(-x_qdp + x_qp)*std::pow(x_qp - xl, -2.0) - 1)/T_q0p;
  J[1*ndim + 3] = h*(x_q - x_qp)*(-x_qdp + x_qp)*std::pow(x_qp - xl, -2.0)/T_q0p;
  J[1*ndim + 8] = -h*(x_q - x_qp)*(-(-x_qdp + x_qp)*std::pow(x_qp - xl, -1.0) + 1)/T_q0p;

  J[2*ndim + 0] = -h/T_d0dp;
  J[2*ndim + 2] = 1 + h/T_d0dp;
  J[2*ndim + 9] = -h*(-x_dp + xl)/T_d0dp;

  J[3*ndim + 1] = h/T_q0dp;
  J[3*ndim + 3] = 1 + h/T_q0dp;
  J[3*ndim + 8] = -h*(-x_qp + xl)/T_q0dp;

  J[4*ndim + 0] = 0.5*h*i_q*(x_ddp - xl)/(H*(x_dp - xl));
  J[4*ndim + 1] = -0.5*h*i_d*(-x_ddp + xl)/(H*(x_qp - xl));
  J[4*ndim + 2] = 0.5*h*i_q*(-x_ddp + x_dp)/(H*(x_dp - xl));
  J[4*ndim + 3] = -0.5*h*i_d*(-x_ddp + x_qp)/(H*(x_qp - xl));
  J[4*ndim + 4] = 1;
  J[4*ndim + 8] = -0.5*h*(-e_qp*(x_ddp - xl)/(x_dp - xl) - phi_1d*(-x_ddp + x_dp)/(x_dp - xl))/H;
  J[4*ndim + 9] = -0.5*h*(e_dp*(-x_ddp + xl)/(x_qp - xl) + phi_2q*(-x_ddp + x_qp)/(x_qp - xl))/H;

  J[5*ndim + 4] = -120.0*M_PI*h;
  J[5*ndim + 5] = 1.0;

  J[6*ndim + 0] = -(x_ddp - xl)/(x_ddp*(x_dp - xl));
  J[6*ndim + 2] = -(-x_ddp + x_dp)/(x_ddp*(x_dp - xl));
  J[6*ndim + 6] = 1/x_ddp;
  J[6*ndim + 9] = 1.0;

  J[7*ndim + 1] = -(-x_qdp + xl)/(x_qdp*(x_qp - xl));
  J[7*ndim + 3] = -(-x_qdp + x_qp)/(x_qdp*(x_qp - xl));
  J[7*ndim + 7] = -1/x_qdp;
  J[7*ndim + 8] = 1.0;

  J[8*ndim + 5] = -v1m*cos(delta - v1a);
  J[8*ndim + 7] = 1.0;
  J[8*ndim + 10] = -sin(delta - v1a);
  J[8*ndim + 11] = v1m*cos(delta - v1a);

  J[9*ndim + 5] = v1m*sin(delta - v1a);
  J[9*ndim + 6] = 1.0;
  J[9*ndim + 10] = -cos(delta - v1a);
  J[9*ndim + 11] = -v1m*sin(delta - v1a);

  J[10*ndim + 6] = i_q;
  J[10*ndim + 7] = i_d;
  J[10*ndim + 8] = v_q;
  J[10*ndim + 9] = v_d;
  J[10*ndim + 10] = v0m*sin(v0a - v1a)/xline;
  J[10*ndim + 11] = -v0m*v1m*cos(v0a - v1a)/xline;

  J[11*ndim + 6] = i_d;
  J[11*ndim + 7] = -i_q;
  J[11*ndim + 8] = -v_d;
  J[11*ndim + 9] = v_q;
  J[11*ndim + 10] = v0m*cos(v0a - v1a)/xline - 2.0*v1m/xline;
  J[11*ndim + 11] = v0m*v1m*sin(v0a - v1a)/xline;

}

void integrate(active* x, size_t dim, System* sys, double h) {
  
  double eps = 1e-9;
  int iteration = 0;
  active *xold = new active[dim];
  active *y = new active[dim];
  for (size_t i = 0; i < dim; ++i) {
    xold[i] = x[i];
  }
  residual_beuler(x, xold, sys, h, y);
  
  active **J=new active*[dim];
  J[0] = new active[dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    J[i] = J[0]+dim*i;
  }
  
  jac_beuler<active>(x, xold, sys, 0.0004, &J[0][0]);
  // linearize Jacobian for BLAS
  double **blasJ=new double*[dim];
  blasJ[0] = new double[dim*dim];
  std::cout << "J going into BLAS" << std::endl;
  for(size_t i = 0; i < dim; ++i) {
    blasJ[i] = blasJ[0]+dim*i;
    for(size_t j = 0; j < dim; ++j) {
      blasJ[i][j]=J[i][j].getValue();
      std::cout << blasJ[i][j] << " ";
    }
    std::cout << std::endl;
  }
  // linearize the RHS for BLAS
  double *py=new double[dim];
  for(size_t i = 0; i < dim; ++i) {
    py[i] = y[i].getValue();
  }
  std::cout << std::endl;
  int ierr=solve(blasJ,py,dim);
  if(ierr) {
    std::cout << "Linear solver error: " << ierr << std::endl;
    exit(1);
  }
  for(size_t i = 0; i < dim; ++i) {
    y[i] = py[i];
  }
  std::cout << "New x and step" << std::endl;
  for(size_t i = 0; i < dim; ++i) {
    x[i]=xold[i]-y[i];
    std::cout << x[i] << " " << y[i] << " " << std::endl;
  }
  std::cout << std::endl;
}

void jactest(double* xold, int dim, System* sys, int h) {
  active *x = new active[dim];
  active *axold = new active[dim];
  active *y = new active[dim];
  for (size_t i = 0; i < dim; ++i) {
    axold[i] = xold[i];
    x[i] = xold[i];
  }

  // Evaluate jacobian
  double J[dim][dim];

  for (size_t j = 0; j < dim; ++j) {
    x[j].setGradient(1.0);
    for (size_t i = 0; i < dim; ++i) {
      residual_beuler(x, axold, sys, h, y);
      J[i][j] = y[i].getGradient();
      //std::cout <<"df/dx: " << y[i].getGradient() << std::endl;
    }
    x[j].setGradient(0.0);
  }

  std::cout << "AD Jacobian" << std::endl;
  // Print jacobian
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      std::cout << J[i][j] << " ";
    }
    std::cout << std::endl;
  }

  std::cout << "HC Jacobian" << std::endl;
  // Hand coded jacobian
  active Jhc[12*12];
  jac_beuler(x, axold, sys, h, Jhc);
  // Print jacobian
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      std::cout << Jhc[i*12 + j] << " ";
    }
    std::cout << std::endl;
  }
  
}

int main(int nargs, char** args) {
  
  
  // Define state arrays
  size_t dim = 12;
  double h = 0.0004;
  
  active *axold=new active[dim];
  double *xold=new double[dim];
  
  
  // System parameters.

  System sys;

  sys.x_d = 1.575; 
  sys.x_q = 1.512;
  sys.x_dp = 0.29;
  sys.x_qp = 0.39 ;
  sys.x_ddp = 0.1733; 
  sys.x_qdp = 0.1733;
  sys.xl = 0.0787;
  sys.H = 3.38; 
  sys.T_d0p = 6.09;
  sys.T_q0p = 1.0;
  sys.T_d0dp = 0.05;
  sys.T_q0dp = 0.15;
  sys.v0m = 1.0162384;
  sys.v0a = -0.05807256;
  sys.xline = 0.0576;


  // Initial values for state array.

  xold[0] = 1.06512;
  xold[1] = 0.5182;
  xold[2] = 0.850584;
  xold[3] = -0.661975;
  xold[4] = -0.05;
  xold[5] = 0.736183;
  xold[6] = 0.770678;
  xold[7] = 0.698323;
  xold[8] = 0.461854;
  xold[9] = 1.01532;
  xold[10] = 1.04;
  xold[11] = 0.0;
  
  jactest(xold, dim, &sys, h);
  for(int i = 0; i < dim ; ++i) axold[i] = xold[i];
  integrate(axold, dim, &sys, h);
  
  return 0;
}
