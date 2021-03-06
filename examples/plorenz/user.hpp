/*!
   \file "user.hpp"
   \brief "This is the provided user code."
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/
#ifndef USER_HPP
#define USER_HPP
#include <iostream>
#include <alg.hpp>

using namespace std;

using namespace alg;

typedef struct System {
public:
  size_t dimension = 7; /* Number of variables */
  double F = 2.0; /* Forcing constant */
  const double h = 0.001;
  
void ic(pVector<double> &x) {
  // Initial values for state array.
  for(size_t i = 0; i < dimension ; ++i) x[i] = F;
  // initial perturbation
  x[0] += 0.01;
}

size_t dim() { return dimension; }

template <class T> void residual_beuler(const pVector<T> &x, const pVector<T> &xold,
    pVector<T> &F) {

  const size_t N = dim(); 
  const double Force = this->F;

  F[0] = (x[1] - x[N - 2])*x[N - 1] - x[0];
  F[1] = (x[2] - x[N - 1])*x[0] - x[1];
  F[N-1] = (x[0] - x[N-3])*x[N - 2] - x[N - 1];

  for (size_t  i = 0; i < N - 3; ++i)
      F[i + 2] = (x[i + 3] - x[i])*x[i + 1] - x[i + 2];

  for (size_t  i = 0; i < N; ++i) {
    F[i] += Force;
    F[i] = x[i] - xold[i] - h*F[i];
  }

}

template <class T> void jac_beuler(const pVector<T> &x, 
    const pVector<T> &xold, pMatrix<T> &J) {
    
  const size_t N = dim(); 

  // REMEMBER the jacobian indeces are switched
  

  // Fill diagonals with ones
  for (size_t  i = 0; i < N; ++i)
    J[i][i] = 1.0;

  // Add integration term
  for (size_t  i = 0; i < N; ++i) {
    if (i == 0) {
      J[i][i] += h;
      J[i][N - 2] -= -h*x[N - 1];
      J[i][N - 1] -= h*(x[i + 1] - x[N - 2]);
      J[i][i + 1] -= h*x[N - 1];
    } else if (i == 1) {
      J[i][i] += h;
      J[i][N - 1] -= -h*x[0];
      J[i][0] -= h*(x[i + 1] - x[N - 1]);
      J[i][i + 1] -= h*x[0];
    } else if (i == (N - 1)) {
      J[i][i] += h;
      J[i][i - 2] -= -h*x[i - 1];
      J[i][i - 1] -= h*(x[0] - x[i - 2]);
      J[i][0] -= h*x[i - 1];
    } else {
      J[i][i] += h;
      J[i][i - 2] -= -h*x[i - 1];
      J[i][i - 1] -= h*(x[i + 1] - x[i - 2]);
      J[i][i + 1] -= h*x[i - 1];
    }
  }

}
} System;
#endif
