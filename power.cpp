#include <math.h>
#include <codi.hpp>
#include <iostream>

void residual_beuler(const codi::RealForward* x, double* xold,
    codi::RealForward* F) {

  // (TEMP): Just put all the parameters here for now.

  double h = 0.0;
  double e_fd = 0.0;
  double p_m = 0.0;

  double x_d = 1.575; 
  double x_q = 1.512;
  double x_dp = 0.29;
  double x_qp = 0.39 ;
  double x_ddp = 0.1733; 
  double x_qdp = 0.1733;
  double xl = 0.0787;
  double H = 3.38; 
  double T_d0p = 6.09;
  double T_q0p = 1.0;
  double T_d0dp = 0.05;
  double T_q0dp = 0.15;

  // infinite bus
  double v0m = 1.0162384;
  double v0a = -0.05807256;
  double xline = 0.0576;

  codi::RealForward e_qp     = x[1];
  codi::RealForward e_dp     = x[2];
  codi::RealForward phi_1d   = x[3];
  codi::RealForward phi_2q   = x[4];
  codi::RealForward w        = x[5];
  codi::RealForward delta    = x[6];
  codi::RealForward v_q      = x[7];
  codi::RealForward v_d      = x[8];
  codi::RealForward i_q      = x[9];
  codi::RealForward i_d      = x[10];
  codi::RealForward v1m      = x[11];
  codi::RealForward v1a      = x[12];

  // Auxiliary variables.
  codi::RealForward psi_de, psi_qe;

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


int main(int nargs, char** args) {
  
  size_t dim = 12;
  
  codi::RealForward x[dim];
  codi::RealForward y[dim];
  double xold[dim];
  
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
  
  for (size_t i = 0; i < dim; ++i)
    x[i] = xold[i];

  // Evaluate something

  x[7].setGradient(1.0);
  for (size_t i = 0; i < dim; ++i) {
    residual_beuler(x, xold, y);
    std::cout <<"df/dx: " << y[i].getGradient() << std::endl;
  }
  x[7].setGradient(0.0);

  return 0;
}
