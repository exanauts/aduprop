#include <cmath>
#include <codi.hpp>
#include <iostream>
#include "linsolve.hpp"

using namespace std;

typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

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

void t1s_driver(double* xic, size_t dim, System* sys, double h, double** J);

template <class T> void residual_beuler(const T* const x, const T* const xold,
    System* const sys, const double h, T* const F) {

  // (TEMP): Just put all the parameters here for now.

  const double e_fd = 2.36980307364616349375;
  const double p_m = 1.06496000000000012875;

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


template <class T> void jac_beuler(const T* const x, const T* const xold,
    System* const sys, const double h, T** const J) {
    

  size_t ndim = 12;

  const double e_fd = 2.36980307364616349375;
  const double p_m = 1.06496000000000012875;

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

void t2s_t1s_integrate(t2s* x, size_t dim, System* sys, double h) {
  
  double eps = 1e-9;
  int iteration = 0;
  t2s *xold = new t2s [dim];
  t2s *y = new t2s [dim];
  for (size_t i = 0; i < dim; ++i) {
    xold[i] = x[i];
  }
  
  residual_beuler<t2s>(x, xold, sys, h, y);
  
  t2s **J = new t2s* [dim];
  J[0] = new t2s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    J[i] = J[0] + dim * i;
  }
  for(size_t i = 0 ; i < dim*dim ; ++i) J[0][i]=0;
  
  jac_beuler<t2s>(x, xold, sys, h, J);
  
  // Get the values and tangents out for both J and y
  double *py = new double[dim];
  double *t1_py = new double[dim];
  double *t2_py = new double[dim];
  double *t2_t1_py = new double[dim];
  double **pJ = new double* [dim];
  double **t1_pJ = new double* [dim];
  double **t2_pJ = new double* [dim];
  double **t2_t1_pJ = new double* [dim];
  pJ[0] = new double [dim*dim];
  t1_pJ[0] = new double [dim*dim];
  t2_pJ[0] = new double [dim*dim];
  t2_t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    pJ[i] = pJ[0] + dim * i;
    t1_pJ[i] = t1_pJ[0] + dim * i;
    t2_pJ[i] = t2_pJ[0] + dim * i;
    t2_t1_pJ[i] = t2_t1_pJ[0] + dim * i;
    py[i] = y[i].value().value();
    t1_py[i] = y[i].value().gradient();
    t2_py[i] = y[i].gradient().value();
    t2_t1_py[i] = y[i].gradient().gradient();
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].value().value();
      t1_pJ[i][j] = J[i][j].value().gradient();
      t2_pJ[i][j] = J[i][j].gradient().value();
      t2_t1_pJ[i][j] = J[i][j].gradient().gradient();
    }
  }

  // Solve 1st order system for t1 and t2
  int ierr = solve(pJ,py,dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x 
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already 
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  decmatmul(t2_pJ, py, t2_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].value().value();
    }
  }
  ierr = solve(pJ,t1_py,dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].value().value();
    }
  }
  ierr = solve(pJ,t2_py,dim);
  decmatmul(t2_t1_pJ, py, t2_t1_py, dim);
  decmatmul(t1_pJ, t2_py, t2_t1_py, dim);
  decmatmul(t2_pJ, t1_py, t2_t1_py, dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].value().value();
    }
  }
  ierr = solve(pJ,t2_t1_py,dim);
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for(size_t i = 0; i < dim; ++i) {
    y[i] = py[i];
    y[i].gradient().gradient() = t2_t1_py[i];
    y[i].value().gradient() = t1_py[i];
    y[i].gradient().value() = t2_py[i];
    x[i]=x[i]-y[i];
    // cout << "py: " << py[i] << " t2_t1_py: " << t2_t1_py[i];
    // cout << " t1_py: " << t1_py[i] << " t2_py: " << t2_py[i];
    // cout << endl;
  }
  // cout << endl;
  delete [] xold; delete [] y; delete [] py; delete [] t1_py;
  delete [] t2_py; delete [] t2_t1_py; delete [] pJ[0]; delete [] t1_pJ[0];
  delete [] t2_pJ[0]; delete [] t2_t1_pJ[0]; delete [] pJ; delete [] t1_pJ;
  delete [] t2_pJ; delete [] t2_t1_pJ; delete [] J[0]; delete [] J;
}

void t1s_integrate(t1s* x, size_t dim, System* sys, double h) {
  
  double eps = 1e-9;
  int iteration = 0;
  t1s *xold = new t1s [dim];
  t1s *y = new t1s [dim];
  for (size_t i = 0; i < dim; ++i) {
    xold[i] = x[i];
  }
  
  residual_beuler<t1s>(x, xold, sys, h, y);
  
  t1s **J = new t1s* [dim];
  J[0] = new t1s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    J[i] = J[0] + dim * i;
  }
  for(size_t i = 0 ; i < dim*dim ; ++i) J[0][i]=0;
  
  jac_beuler<t1s>(x, xold, sys, h, J);
  
  // Get the values and tangents out for both J and y
  double *py = new double[dim];
  double *t1_py = new double[dim];
  double **pJ = new double* [dim];
  double **t1_pJ = new double* [dim];
  pJ[0] = new double [dim*dim];
  t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    pJ[i] = pJ[0] + dim * i;
    t1_pJ[i] = t1_pJ[0] + dim * i;
    py[i] = y[i].getValue();
    t1_py[i] = y[i].getGradient();
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].getValue();
      t1_pJ[i][j] = J[i][j].getGradient();
    }
  }

  // Solve 1st order system
  int ierr = solve(pJ,py,dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x 
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already 
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = J[i][j].getValue();
    }
  }
  ierr = solve(pJ,t1_py,dim);
  // That's it, we have the tangents t1_x in t1_py of the LS Ax=b
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for(size_t i = 0; i < dim; ++i) {
    y[i] = py[i];
    y[i].setGradient(t1_py[i]);
    x[i]=x[i]-y[i];
    // cout << x[i] << " " << y[i] << " " << endl;
  }
  // cout << endl;
  delete [] xold; delete [] y; delete [] py; delete [] t1_py;
  delete [] pJ[0]; delete [] pJ; delete [] t1_pJ[0];
  delete [] t1_pJ; delete [] J[0]; delete [] J;
}

void integrate(double* x, size_t dim, System* sys, double h) {
  
  double eps = 1e-9;
  int iteration = 0, ierr;
  double *xold = new double[dim];
  double *y = new double[dim];
  double **J = new double* [dim];
  
  for (size_t i = 0; i < dim; ++i)
    xold[i] = x[i];
  
  residual_beuler<double>(x, xold, sys, h, y);
  
  J[0] = new double [dim*dim];
  
  for(size_t i = 0; i < dim; ++i)
    J[i] = J[0] + dim * i;

  for(int i = 0; i < dim*dim; ++i)
    J[0][i] = 0;
  
  jac_beuler<double>(x, xold, sys, h, J);
  
  ierr = solve(J, y, dim);
  
  if(ierr) {
    cout << "Linear solver error: " << ierr << endl;
    exit(1);
  }

  //cout << "New x and step" << endl;
  for(size_t i = 0; i < dim; ++i) {
    x[i] = x[i] - y[i];
    //cout << x[i] << " = " << xold[i] << " - " << y[i] << " " << endl;
  }

  // cout << endl;
  delete [] xold;
  delete [] y;
  delete [] J[0];
  delete [] J;
}

void jactest(double* xold, size_t dim, System* sys, double h) {
  t1s *x = new t1s[dim];
  t1s *axold = new t1s[dim];
  t1s *y = new t1s[dim];
  for (size_t i = 0; i < dim; ++i) {
    axold[i] = xold[i];
    x[i] = xold[i];
  }

  // Evaluate jacobian
  double** J = new double*[dim];
  J[0] = new double[dim*dim];
  for(size_t i = 0; i < dim ; ++i) J[i] = J[0] + i * dim;

  for (size_t j = 0; j < dim; ++j) {
    x[j].setGradient(1.0);
    for (size_t i = 0; i < dim; ++i) {
      residual_beuler(x, axold, sys, h, y);
      J[i][j] = y[i].getGradient();
    }
    x[j].setGradient(0.0);
  }

  cout << "AD Jacobian" << endl;
  // Print jacobian
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      cout << J[i][j] << " ";
    }
    cout << endl;
  }

  cout << "HC Jacobian" << endl;
  // Hand coded jacobian
  t1s **Jhc= new t1s*[dim];
  Jhc[0] = new t1s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    Jhc[i] = Jhc[0] + dim * i;
  }
  for(size_t i = 0 ; i < dim*dim ; ++i) Jhc[0][i]=0;
  
  
  jac_beuler(x, axold, sys, h, Jhc);
  // Print jacobian
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      cout << Jhc[i][j] << " ";
    }
    cout << endl;
  }
  delete [] x;
  delete [] axold;
  delete [] y;
}

// Driver for accumulating Jacobian using FD
void fdJ_driver(double* xic, size_t dim, System* sys, double h, double** J) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double pert=1e-8;
  
  for(size_t i = 0 ; i < dim ; ++i) {
    for(size_t j = 0 ; j < dim ; ++j) xpert1[j] = xic[j];
    for(size_t j = 0 ; j < dim ; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    integrate(xpert1, dim, sys, h);
    integrate(xpert2, dim, sys, h);
    for(size_t j = 0 ; j < dim ; ++j) J[i][j]=(xpert1[j]-xpert2[j])/pert;
  }
  integrate(xic, dim, sys, h);
  delete [] xpert1;  
  delete [] xpert2;  
}

void fdH_driver(double* xic, size_t dim, System* sys, double h, double*** H) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double **Jpert1 = new double*[dim];
  Jpert1[0] = new double[dim*dim];
  for(size_t i = 0; i < dim; ++i) Jpert1[i] = Jpert1[0] + i*dim;
  double **Jpert2 = new double*[dim];
  Jpert2[0] = new double[dim*dim];
  for(size_t i = 0; i < dim; ++i) Jpert2[i] = Jpert2[0] + i*dim;
  double pert=1e-8;
  
  for(size_t i = 0 ; i < dim ; ++i) {
    for(size_t j = 0 ; j < dim ; ++j) xpert1[j] = xic[j];
    for(size_t j = 0 ; j < dim ; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t1s_driver(xpert1, dim, sys, h, Jpert1);
    t1s_driver(xpert2, dim, sys, h, Jpert2);
    for(size_t j = 0 ; j < dim ; ++j) {
      for(size_t k = 0 ; k < dim ; ++k) {
        H[i][j][k] = (Jpert1[j][k]-Jpert2[j][k])/pert;
      } 
    }
  }
  delete [] xpert1;    delete [] xpert2;
  delete [] Jpert1[0]; delete [] Jpert1;
  delete [] Jpert2[0]; delete [] Jpert2;
}

// Driver for accumulating Jacobian using AD
void t1s_driver(double* xic, size_t dim, System* sys, double h, double** J) {
  t1s* axic = new t1s [dim];
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim ; ++j) {
      axic[j] = xic[j];
      axic[j].setGradient(0.0);
    }
    axic[i].setGradient(1.0);
    t1s_integrate(axic, dim, sys, h);
    for(size_t j = 0; j < dim ; ++j) J[i][j] = axic[j].getGradient();
  }
  for(size_t i = 0; i < dim ; i++) xic[i] = axic[i].getValue();
  delete [] axic;
}

void t2s_t1s_driver(double* xic, size_t dim, System* sys, int h, double **J, double*** H) {
  // t2s* axic = new t2s [dim];
  for(size_t i = 0; i < dim; ++i) {
    t2s* axic = new t2s [dim];
    for(size_t j = 0; j < dim ; ++j) {
      for(size_t k = 0; k < dim ; ++k) {
        axic[k].value().value()       = xic[k];
        axic[k].gradient().gradient() = 0.0;
        axic[k].value().gradient()    = 0.0;
        axic[k].gradient().value()    = 0.0;
      }
      axic[i].value().gradient() = 1.0;
      axic[j].gradient().value() = 1.0;
      t2s_t1s_integrate(axic, dim, sys, h);
      for(size_t k = 0; k < dim ; ++k) {
        H[i][j][k] = axic[k].gradient().gradient();
      }
      // This should give you the Jacobian too
      // for(size_t k = 0; k < dim ; ++k) {
      //   J[j][k] = axic[k].gradient().value();
      // }
    }
    for(size_t j = 0; j < dim ; ++j) {
      J[i][j] = axic[j].value().gradient();
    }
    delete [] axic;
  }
  // for(size_t i = 0; i < dim ; i++) xic[i] = axic[i].value().value();
  // delete [] axic;
}

int main(int nargs, char** args) {
  
  
  // Define state arrays
  size_t dim = 12;
  double h = 0.004;
  
  double *xold = new double [dim];
  double *x = new double [dim];
  
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

  xold[0] = 1.06512037300928485983;
  xold[1] = 0.51819992367912581788;
  xold[2] = 0.85058383242985102779;
  xold[3] = -0.66197500054304025952;
  xold[4] = -0.05000000000000000278;
  xold[5] = 0.73618306350367335167;
  xold[6] = 0.77067836274882195458;
  xold[7] = 0.69832289180288620312;
  xold[8] = 0.46185376441989828278;
  xold[9] = 1.01531727676021699125;
  xold[10] = 1.0400000000000000000;
  xold[11] = 0.00000000000000000000;
  
  // jactest(xold, dim, &sys, h);
  cout << "At point:" << endl;
  cout << "---------" << endl;
  for(size_t i = 0; i < dim; ++i) cout << xold[i] << " ";
  cout << endl;
  
  double** J = new double*[dim]; 
  J[0] = new double[dim*dim];
  for(size_t i = 0; i < dim; ++i) J[i] = J[0] + i * dim;
  double ***H = new double** [dim];
  cout << "Total: " << dim*dim*dim << endl;
  H[0] = new double* [dim*dim];
  H[0][0] = new double [dim*dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    H[i] = H[0] + i * dim;
    for(size_t j = 0; j < dim; ++j) {
      H[i][j] = H[0][0] + i * dim * dim + j * dim;
    }
  }
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  t1s_driver(x, dim, &sys, h, J);
  cout << "Function using AD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    cout << xold[i] << " ";
  }
  cout << endl;
  cout << endl;
  cout << "Jacobian using 1st order AD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      cout << J[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  cout << "Hessian and Jacobian using 2nd order AD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      J[i][j] = 0.0;
    }  
  }
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  t2s_t1s_driver(x, dim, &sys, h, J, H);
  cout << "J" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      cout << J[i][j] << " ";
    }
    cout << endl;
  }
  cout << "H" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      for(size_t k = 0; k < dim; ++k) {
        cout << H[i][j][k] << " ";
      }
    }
    cout << endl;
  }
  cout << endl;
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  fdJ_driver(x, dim, &sys, h, J);
  cout << "Function using FD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    cout << xold[i] << " ";
  }
  cout << endl;
  cout << endl;
  cout << "Jacobian using FD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      cout << J[i][j] << " ";
    }
    cout << endl;
  }
  
  cout << endl;
  cout << "Hessian using FD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  fdH_driver(x, dim, &sys, h, H);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      for(size_t k = 0; k < dim; ++k) {
        cout << H[i][j][k] << " ";
      }
    }
    cout << endl;
  }
  cout << endl;
  delete [] H[0][0];
  delete [] H[0];
  delete [] H;
  delete [] J[0];
  delete [] J;
  delete [] xold;
  delete [] x;
  return 0;
}
