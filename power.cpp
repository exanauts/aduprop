/*!
   \file "power.cpp"
   \brief "Jacobian, Hessian and Tensor accumulation using Automatic
           Differentiation."
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <cmath>
// AD library used
#include <codi.hpp>
#include <iostream>
// Rudimentary linear solver interface
#include "ad.hpp"

using namespace std;



int main(int nargs, char** args) {
  
  // Define state arrays
  size_t dim = 12;
  init(dim);
  double h = 0.004;
  
  double *xold = new double [dim];
  double *x = new double [dim];
  
  // System parameters.
  
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
  H[0] = new double* [dim*dim];
  H[0][0] = new double [dim*dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    H[i] = H[0] + i * dim;
    for(size_t j = 0; j < dim; ++j) {
      H[i][j] = H[0][0] + i * dim * dim + j * dim;
    }
  }
  double ****T = new double*** [dim];
  T[0] = new double** [dim*dim];
  T[0][0] = new double* [dim*dim*dim];
  T[0][0][0] = new double [dim*dim*dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    T[i] = T[0] + i * dim;
    for(size_t j = 0; j < dim; ++j) {
      T[i][j] = T[0][0] + i * dim * dim + j * dim;
      for(size_t k = 0; k < dim; ++k) {
        T[i][j][k] = T[0][0][0] + i * dim * dim * dim + j * dim * dim + k * dim;
      }
    }
  }
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  t1s_driver(x, dim, h, J);
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
  t2s_t1s_driver(x, dim, h, J, H);
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
  cout << "Tensor, Hessian and Jacobian using 3rd order AD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      J[i][j] = 0.0;
    }  
  }
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  t3s_t2s_t1s_driver(x, dim, h, J, H, T);
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
  cout << "T" << endl;
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      for(size_t k = 0; k < dim; ++k) {
        for(size_t l = 0; l < dim; ++l) {
          cout << T[i][j][k][l] << " ";
        }
      }
    }
    cout << endl;
  }
  cout << endl;
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  fdJ_driver(x, dim, h, J);
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
  fdH_driver(x, dim, h, H);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      for(size_t k = 0; k < dim; ++k) {
        cout << H[i][j][k] << " ";
      }
    }
    cout << endl;
  }
  cout << endl;
  cout << "Tensor using FD" << endl;
  cout << "-----------------" << endl;
  for(size_t i = 0; i < dim; ++i) x[i] = xold[i];
  fdT_driver(x, dim, h, T);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      for(size_t k = 0; k < dim; ++k) {
        for(size_t l = 0; l < dim; ++l) {
          cout << T[i][j][k][l] << " ";
        }
      }
    }
    cout << endl;
  }
  cout << endl;
  destroy();
  delete [] T[0][0][0];
  delete [] T[0][0];
  delete [] T[0];
  delete [] T;
  delete [] H[0][0];
  delete [] H[0];
  delete [] H;
  delete [] J[0];
  delete [] J;
  delete [] xold;
  delete [] x;
  return 0;
}
