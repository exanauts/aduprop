/*!
   \file "power.cpp"
   \brief "Jacobian, Hessian and Tensor accumulation using Automatic
           Differentiation. This is the main function where everything comes 
           together"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <cmath>
#include <codi.hpp>
#include <iostream>
// AD library used
// Rudimentary linear solver interface
#include "ad.hpp"

using namespace std;



int main(int nargs, char** args) {
  // Define state arrays
  size_t dim = 12;
  init(dim);
  double h = 0.004;

  pVector<double> xold(dim);
  pVector<double> x(dim);

  // System parameters.
  sys.x_d = 1.575;
  sys.x_q = 1.512;
  sys.x_dp = 0.29;
  sys.x_qp = 0.39;
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

  jactest(xold, dim, h);
  cout << "At point:" << endl;
  cout << "---------" << endl;
  cout << xold << endl;
  pMatrix<double>  J(dim,dim);
  pTensor3<double> H(dim,dim,dim);
  pTensor4<double> T(dim,dim,dim,dim);
  x = xold;
  t1s_driver(x, dim, h, J);
  cout << "Function using AD" << endl;
  cout << "-----------------" << endl;
  cout << xold << endl;
  cout << endl;
  cout << "Jacobian using 1st order AD" << endl;
  cout << "-----------------" << endl;
  cout << J << endl;
  cout << "Hessian and Jacobian using 2nd order AD" << endl;
  cout << "-----------------" << endl;
  
  J.zeros();
  H.zeros();
  x = xold;
  t2s_t1s_driver(x, dim, h, J, H);
  cout << "J" << endl << J << endl;
  cout << "H" << endl << H << endl;
  cout << "Tensor, Hessian and Jacobian using 3rd order AD" << endl;
  cout << "-----------------" << endl;
  J.zeros();
  H.zeros();
  x = xold;
  t3s_t2s_t1s_driver(x, dim, h, J, H, T);
  cout << "J" << endl << J << endl;
  cout << "H" << endl << H << endl;
  cout << "T" << endl << T << endl;
  for (size_t i = 0; i < dim; ++i) x[i] = xold[i];
  fdJ_driver(x, dim, h, J);
  cout << "Function using FD" << endl;
  cout << "-----------------" << endl;
  cout << xold << endl;
  cout << endl;
  cout << "Jacobian using FD" << endl;
  cout << "-----------------" << endl;
  cout << "J" << endl << J << endl;
  cout << "Hessian using FD" << endl;
  cout << "-----------------" << endl;
  x = xold;
  fdH_driver(x, dim, h, H);
  cout << H << endl;
  cout << "Tensor using FD" << endl;
  cout << "-----------------" << endl;
  x = xold;
  fdT_driver(x, dim, h, T);
  cout << T << endl;
  return 0;
}
