#ifndef USER_HPP
#define USER_HPP
#include <iostream>

using namespace std;

template <class T> int adlinsolve(T **A, T *B, size_t n);

typedef struct System {
  size_t N; /* Number of variables */
  double F; /* Forcing constant */
};

System sys;

template <class T> void residual_beuler(const T* const x, const T* const xold,
    const double h, T* const F) {

  const size_t N = sys.N; 
  const double Force = sys.F;

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


template <class T> void jac_beuler(const T* const x, const T* const xold,
    const double h, T** const J) {
    

  const size_t N = sys.N; 
  const double Force = sys.F;

  // REMEMBER the jacobian indeces are switched
  

  // Fill diagonals with ones
  for (size_t  i = 0; i < N; ++i)
    J[i][i] = 1.0;

  // Add integration term
  for (size_t  i = 0; i < N; ++i) {
    if (i == 0) {
      J[i][i] -= h;
      J[N - 2][i] += -h*x[N - 1];
      J[N - 1][i] += h*(x[i + 1] - x[N - 2]);
      J[i + 1][i] += h*x[N - 1];
    } else if (i == 1) {
      J[i][i] -= h;
      J[N - 1][i] += -h*x[0];
      J[0][i] += h*(x[i + 1] - x[N - 1]);
      J[i + 1][i] += h*x[0];
    } else if (i == (N - 1)) {
      J[i][i] -= h;
      J[i - 2][i] += -h*x[i - 1];
      J[i - 1][i] += h*(x[0] - x[i - 2]);
      J[0][i] += h*x[i - 1];
    } else {
      J[i][i] -= h;
      J[i - 2][i] += -h*x[i - 1];
      J[i - 1][i] += h*(x[i + 1] - x[i - 2]);
      J[i + 1][i] += h*x[i - 1];
    }
  }

}

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
