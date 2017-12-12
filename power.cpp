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

  System sys;
  size_t dim = sys.dim();  

  pVector<double> xold(dim);
  pVector<double> x(dim);
  
  sys.ic(xold);
  ad drivers(sys);

  // drivers.jactest(xold);
  cout << "At point:" << endl;
  cout << "---------" << endl;
  cout << xold << endl;
  pMatrix<double>  J(dim,dim);
  pTensor3<double> H(dim,dim,dim);
  pTensor4<double> T(dim,dim,dim,dim);
  x = xold;
  drivers.integrate<double>(x);
  // drivers.t1s_driver(x, J);
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
  // drivers.t2s_t1s_driver(x, J, H);
  // cout << "J" << endl << J << endl;
  // cout << "H" << endl << H << endl;
  cout << "Tensor, Hessian and Jacobian using 3rd order AD" << endl;
  cout << "-----------------" << endl;
  J.zeros();
  H.zeros();
  x = xold;
  // drivers.t3s_t2s_t1s_driver(x, J, H, T);
  // cout << "J" << endl << J << endl;
  // cout << "H" << endl << H << endl;
  // cout << "T" << endl << T << endl;
  for (size_t i = 0; i < dim; ++i) x[i] = xold[i];
  drivers.fdJ_driver(x, J);
  cout << "Function using FD" << endl;
  cout << "-----------------" << endl;
  cout << xold << endl;
  cout << endl;
  cout << "Jacobian using FD" << endl;
  cout << "-----------------" << endl;
  // cout << "J" << endl << J << endl;
  cout << "Hessian using FD" << endl;
  cout << "-----------------" << endl;
  x = xold;
  // drivers.fdH_driver(x, H);
  // cout << H << endl;
  cout << "Tensor using FD" << endl;
  cout << "-----------------" << endl;
  x = xold;
  // drivers.fdT_driver(x, T);
  // cout << T << endl;
  return 0;
}
