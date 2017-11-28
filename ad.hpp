#ifndef LINSOLVE_HPP
#define LINSOLVE_HPP
#include <iostream>
#include "user.hpp"

using namespace std;

typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

void t1s_driver(double* xic, size_t dim, double h, double** J);
void t2s_t1s_driver(double* xic, size_t dim, double h, double **J, double*** H);
void t3s_t2s_t1s_driver(double* xic, size_t dim, double h, double **J, double*** H, double ****T);


double *py; 
double *t3_py;
double *t1_py;
double *t3_t1_py;
double *t2_py;
double *t3_t2_py;
double *t2_t1_py;
double *t3_t2_t1_py;
double **pJ;
double **t3_pJ;
double **t1_pJ;
double **t3_t1_pJ;
double **t2_pJ;
double **t3_t2_pJ;
double **t2_t1_pJ;
double **t3_t2_t1_pJ;
double **tmp_pJ;

t3s *t3s_xold; t3s *t3s_y;
t3s **t3s_J;

t2s *t2s_xold; t2s *t2s_y;
t2s **t2s_J;

t1s *t1s_xold; t1s *t1s_y;
t1s **t1s_J;

#ifndef FNAME
#ifndef __bg__
#define FNAME(f) f ## _ 
#else
#define FNAME(f) f // no underscores for fortran names on bgp
#endif
#endif

// declarations for LAPACK and BLAS functions used to factor/solve:

extern "C" void FNAME(dgesv)(int *n,
      int *nrhs,
      double A[],
      int *lda,
			int ipiv[], 
      double B[],
      int *ldb,
      int *info); // should be 0 on exit
      
extern "C" void FNAME(dsytrf)(char *uplo, 
			int *n, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double work[],
			int *lwork, 
			int *info);

// dsytrs_() solves the system Ax = b using the factor obtained by dsytrf_().
extern "C" void FNAME(dsytrs)(char *uplo, 
			int *n, 
			int *nrhs, 
			double A[], 
			int *lda, 
			int ipiv[], 
			double b[], 
			int *ldb,
			int *info);
extern "C" void FNAME(dgemv)(char *trans, 
			int *m, 
			int *n, 
      double *alpha,
			double A[], 
			int *lda, 
      double *x,
      int *incx,
      double *beta,
      double *y,
      int *incy);
      
void init(size_t dim) {
  py = new double[dim];
  t3_py = new double[dim];
  t1_py = new double[dim];
  t3_t1_py = new double[dim];
  t2_py = new double[dim];
  t3_t2_py = new double[dim];
  t2_t1_py = new double[dim];
  t3_t2_t1_py = new double[dim];
  pJ = new double* [dim];
  t3_pJ = new double* [dim];
  t1_pJ = new double* [dim];
  t3_t1_pJ = new double* [dim];
  t2_pJ = new double* [dim];
  t3_t2_pJ = new double* [dim];
  t2_t1_pJ = new double* [dim];
  t3_t2_t1_pJ = new double* [dim];
  pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    pJ[i] = pJ[0] + dim * i;
  }
  t3_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t3_pJ[i] = t3_pJ[0] + dim * i;
  }
  t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t1_pJ[i] = t1_pJ[0] + dim * i;
  }
  t3_t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t3_t1_pJ[i] = t3_t1_pJ[0] + dim * i;
  }
  t2_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t2_pJ[i] = t2_pJ[0] + dim * i;
  }
  t3_t2_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t3_t2_pJ[i] = t3_t2_pJ[0] + dim * i;
  }
  t2_t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t2_t1_pJ[i] = t2_t1_pJ[0] + dim * i;
  }
  t3_t2_t1_pJ[0] = new double [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t3_t2_t1_pJ[i] = t3_t2_t1_pJ[0] + dim * i;
  }
  
  t3s_xold = new t3s [dim];
  t3s_y = new t3s [dim];
  t3s_J = new t3s* [dim];
  t3s_J[0] = new t3s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t3s_J[i] = t3s_J[0] + dim * i;
  }
  t2s_xold = new t2s [dim];
  t2s_y = new t2s [dim];
  t2s_J = new t2s* [dim];
  t2s_J[0] = new t2s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t2s_J[i] = t2s_J[0] + dim * i;
  }
  t1s_xold = new t1s [dim];
  t1s_y = new t1s [dim];
  t1s_J = new t1s* [dim];
  t1s_J[0] = new t1s [dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    t1s_J[i] = t1s_J[0] + dim * i;
  }
  tmp_pJ = new double*[dim];
  tmp_pJ[0] = new double[dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    tmp_pJ[i] = tmp_pJ[0] + dim * i;
  }
}

void destroy() {
  delete [] py; delete [] t3_py, delete [] t1_py;
  delete [] t3_t1_py;
  delete [] t2_py; delete [] t2_t1_py; delete [] pJ[0]; delete [] t1_pJ[0];
  delete [] t3_t2_py; delete [] t3_t2_t1_py; delete [] t3_pJ[0]; delete [] t3_t1_pJ[0];
  delete [] t2_pJ[0]; delete [] t2_t1_pJ[0]; delete [] pJ; delete [] t1_pJ;
  delete [] t3_t2_pJ[0]; delete [] t3_t2_t1_pJ[0]; delete [] t3_pJ; delete [] t3_t1_pJ;
  delete [] t2_pJ; delete [] t2_t1_pJ; 
  delete [] t3_t2_pJ; delete [] t3_t2_t1_pJ; 
  
  delete [] tmp_pJ[0];
  delete [] tmp_pJ;
  
  delete [] t3s_xold; 
  delete [] t3s_y;
  delete [] t3s_J[0]; delete [] t3s_J;
  
  delete [] t2s_xold; 
  delete [] t2s_y;
  delete [] t2s_J[0]; delete [] t2s_J;
  
  delete [] t1s_xold; 
  delete [] t1s_y;
  delete [] t1s_J[0]; delete [] t1s_J;
}
// Decremental matmul
void decmatmul(double **A, double *x, double *y, size_t n) {
  // for(size_t i = 0; i < n; ++i) {
  //   for(size_t j = 0; j < n; ++j) {
  //     y[i] -= A[j][i] * x[j];
  //   }
  // }
  int n_=(int) n;
  char trans = 'N';
  double alpha = -1.0;
  double beta = 1.0;
  int incx=1;
  int incy=1;
  FNAME(dgemv)(&trans, &n_, &n_, &alpha, &A[0][0], &n_, x, &incx, &beta, y, &incy);
}
      
int solve(double **A, double *B, size_t n) {
  int n_=(int) n;
  int lda=n_;
  int ldb=n_;
  int nrhs=1;
  int info=0;
  static bool first = true;
  static int *ipiv;
  
  if(first) {
    ipiv = new int[n];
    first = false;
  }
  
  FNAME(dgesv)(&n_, &nrhs, &A[0][0], &lda, ipiv, &B[0], &ldb, &info);
  if(info != 0) {
    cout << "Error in LS. Error code: " << info << endl;
    exit(1);
    
  }
  return info;
}

template <class T> int adlinsolve(T **A, T *B, size_t n) {
  T t;
  cout << "adlinsolve not implemented for this type " << endl;
}

template <> int adlinsolve<double>(double **J, double *y, size_t dim) {
  return solve(J, y, dim);
}

template <> int adlinsolve<t1s>(t1s **t1s_J, t1s *t1s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for(size_t i = 0; i < dim; ++i) {
    py[i] = t1s_y[i].getValue();
    t1_py[i] = t1s_y[i].getGradient();
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t1s_J[i][j].getValue();
      tmp_pJ[i][j] = t1s_J[i][j].getValue();
      t1_pJ[i][j] = t1s_J[i][j].getGradient();
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
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t1_py,dim);
  // That's it, we have the tangents t1_x in t1_py of the LS Ax=b
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for(size_t i = 0; i < dim; ++i) {
    t1s_y[i] = py[i];
    t1s_y[i].setGradient(t1_py[i]);
  }
  return 0;
}

template <> int adlinsolve<t2s>(t2s **t2s_J, t2s *t2s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for(size_t i = 0; i < dim; ++i) {
    py[i] = t2s_y[i].value().value();
    t1_py[i] = t2s_y[i].value().gradient();
    t2_py[i] = t2s_y[i].gradient().value();
    t2_t1_py[i] = t2s_y[i].gradient().gradient();
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t2s_J[i][j].value().value();
      tmp_pJ[i][j] = t2s_J[i][j].value().value();
      t1_pJ[i][j] = t2s_J[i][j].value().gradient();
      t2_pJ[i][j] = t2s_J[i][j].gradient().value();
      t2_t1_pJ[i][j] = t2s_J[i][j].gradient().gradient();
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
      pJ[i][j] = tmp_pJ[i][j]; 
    }
  }
  ierr = solve(pJ,t1_py,dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j]; 
    }
  }
  ierr = solve(pJ,t2_py,dim);
  decmatmul(t2_t1_pJ, py, t2_t1_py, dim);
  decmatmul(t1_pJ, t2_py, t2_t1_py, dim);
  decmatmul(t2_pJ, t1_py, t2_t1_py, dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j]; 
    }
  }
  ierr = solve(pJ,t2_t1_py,dim);
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for(size_t i = 0; i < dim; ++i) {
    t2s_y[i] = py[i];
    t2s_y[i].gradient().gradient() = t2_t1_py[i];
    t2s_y[i].value().gradient() = t1_py[i];
    t2s_y[i].gradient().value() = t2_py[i];
  }
  return 0;
}

template <> int adlinsolve<t3s>(t3s **t3s_J, t3s *t3s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for(size_t i = 0; i < dim; ++i) {
    py[i] = t3s_y[i].value().value().value();
    t3_py[i] = t3s_y[i].gradient().value().value();
    t1_py[i] = t3s_y[i].value().value().gradient();
    t3_t1_py[i] = t3s_y[i].gradient().value().gradient();
    t2_py[i] = t3s_y[i].value().gradient().value();
    t3_t2_py[i] = t3s_y[i].gradient().gradient().value();
    t2_t1_py[i] = t3s_y[i].value().gradient().gradient();
    t3_t2_t1_py[i] = t3s_y[i].gradient().gradient().gradient();
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t3s_J[i][j].value().value().value();
      tmp_pJ[i][j] = t3s_J[i][j].value().value().value();
      t3_pJ[i][j] = t3s_J[i][j].gradient().value().value();
      t1_pJ[i][j] = t3s_J[i][j].value().value().gradient();
      t3_t1_pJ[i][j] = t3s_J[i][j].gradient().value().gradient();
      t2_pJ[i][j] = t3s_J[i][j].value().gradient().value();
      t3_t2_pJ[i][j] = t3s_J[i][j].gradient().gradient().value();
      t2_t1_pJ[i][j] = t3s_J[i][j].value().gradient().gradient();
      t3_t2_t1_pJ[i][j] = t3s_J[i][j].gradient().gradient().gradient();
    }
  }

  // Solve 1st order system for t1 and t2
  int ierr = solve(pJ,py,dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x 
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already 
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  decmatmul(t2_pJ, py, t2_py, dim);
  decmatmul(t3_pJ, py, t3_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t1_py,dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t2_py,dim);
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t3_py,dim);
  
  decmatmul(t2_t1_pJ, py, t2_t1_py, dim);
  decmatmul(t1_pJ, t2_py, t2_t1_py, dim);
  decmatmul(t2_pJ, t1_py, t2_t1_py, dim);
  
  decmatmul(t3_t2_pJ, py, t3_t2_py, dim);
  decmatmul(t3_pJ, t2_py, t3_t2_py, dim);
  decmatmul(t2_pJ, t3_py, t3_t2_py, dim);
  
  decmatmul(t3_t1_pJ, py, t3_t1_py, dim);
  decmatmul(t1_pJ, t3_py, t3_t1_py, dim);
  decmatmul(t3_pJ, t1_py, t3_t1_py, dim);
  
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t2_t1_py,dim);
  
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t3_t2_py,dim);
  
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t3_t1_py,dim);
  
  decmatmul(t3_t1_pJ, t2_py, t3_t2_t1_py, dim);
  decmatmul(t1_pJ, t3_t2_py, t3_t2_t1_py, dim);
  decmatmul(t3_t2_t1_pJ, py, t3_t2_t1_py, dim);
  decmatmul(t2_t1_pJ, t3_py, t3_t2_t1_py, dim);
  decmatmul(t3_t2_pJ, t1_py, t3_t2_t1_py, dim);
  decmatmul(t2_pJ, t3_t1_py, t3_t2_t1_py, dim);
  decmatmul(t3_pJ, t2_t1_py, t3_t2_t1_py, dim);
  
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ,t3_t2_t1_py,dim);
  
  // Put x and t1_x back into the t1s type
  for(size_t i = 0; i < dim; ++i) {
    t3s_y[i].value().value().value() = py[i];
    t3s_y[i].gradient().value().value() = t3_py[i];
    t3s_y[i].value().gradient().gradient() = t2_t1_py[i];
    t3s_y[i].gradient().gradient().gradient() = t3_t2_t1_py[i];
    t3s_y[i].value().value().gradient() = t1_py[i];
    t3s_y[i].gradient().value().gradient() = t3_t1_py[i];
    t3s_y[i].value().gradient().value() = t2_py[i];
    t3s_y[i].gradient().gradient().value() = t3_t2_py[i];
  }
  return 0;
}
// Driver for accumulating Jacobian using FD
void fdJ_driver(double* xic, size_t dim, double h, double** J) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double pert=1e-8;
  
  for(size_t i = 0 ; i < dim ; ++i) {
    for(size_t j = 0 ; j < dim ; ++j) xpert1[j] = xic[j];
    for(size_t j = 0 ; j < dim ; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    integrate<double>(xpert1, dim, h);
    integrate<double>(xpert2, dim, h);
    for(size_t j = 0 ; j < dim ; ++j) J[i][j]=(xpert1[j]-xpert2[j])/pert;
  }
  integrate<double>(xic, dim, h);
  delete [] xpert1;  
  delete [] xpert2;  
}

void fdH_driver(double* xic, size_t dim, double h, double*** H) {
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
    t1s_driver(xpert1, dim, h, Jpert1);
    t1s_driver(xpert2, dim, h, Jpert2);
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

void fdT_driver(double* xic, size_t dim, double h, double**** T) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double **J = new double* [dim];
  J[0] = new double[dim*dim];
  for(size_t i = 0; i < dim; ++i) J[i] = J[0] + i*dim;
  double ***Hpert1 = new double** [dim];
  Hpert1[0] = new double* [dim*dim];
  Hpert1[0][0] = new double [dim*dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    Hpert1[i] = Hpert1[0] + i * dim;
    for(size_t j = 0; j < dim; ++j) {
      Hpert1[i][j] = Hpert1[0][0] + i * dim * dim + j * dim;
    }
  }
  double ***Hpert2 = new double** [dim];
  Hpert2[0] = new double* [dim*dim];
  Hpert2[0][0] = new double [dim*dim*dim];
  for(size_t i = 0; i < dim; ++i) {
    Hpert2[i] = Hpert2[0] + i * dim;
    for(size_t j = 0; j < dim; ++j) {
      Hpert2[i][j] = Hpert2[0][0] + i * dim * dim + j * dim;
    }
  }
  double pert=1e-8;
  for(size_t i = 0 ; i < dim ; ++i) {
    cout << "Computing tensor. " << (double) i*(double) 100.0/(double) dim << "\% done." << endl;
    for(size_t j = 0 ; j < dim ; ++j) xpert1[j] = xic[j];
    for(size_t j = 0 ; j < dim ; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t2s_t1s_driver(xpert1, dim, h, J, Hpert1);
    t2s_t1s_driver(xpert2, dim, h, J, Hpert2);
    for(size_t j = 0 ; j < dim ; ++j) {
      for(size_t k = 0 ; k < dim ; ++k) {
        for(size_t l = 0 ; l < dim ; ++l) {
          T[i][j][k][l] = (Hpert1[j][k][l]-Hpert2[j][k][l])/pert;
        }
      } 
    }
  }
  delete [] Hpert1[0][0];
  delete [] Hpert1[0];
  delete [] Hpert1;
  delete [] Hpert2[0][0];
  delete [] Hpert2[0];
  delete [] Hpert2;
  delete [] J[0];
  delete [] J;
  delete [] xpert1;    delete [] xpert2;
}

// Driver for accumulating Jacobian using AD
void t1s_driver(double* xic, size_t dim, double h, double** J) {
  t1s* axic = new t1s [dim];
  for(size_t i = 0; i < dim; ++i) {
    for(size_t j = 0; j < dim ; ++j) {
      axic[j] = xic[j];
      axic[j].setGradient(0.0);
    }
    axic[i].setGradient(1.0);
    integrate<t1s>(axic, dim, h);
    for(size_t j = 0; j < dim ; ++j) J[i][j] = axic[j].getGradient();
  }
  for(size_t i = 0; i < dim ; i++) xic[i] = axic[i].getValue();
  delete [] axic;
}

void t2s_t1s_driver(double* xic, size_t dim, double h, double **J, double*** H) {
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
      integrate<t2s>(axic, dim, h);
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
}

void t3s_t2s_t1s_driver(double* xic, size_t dim, double h, double **J, double*** H, double ****T) {
  t3s* axic = new t3s [dim];
  for(size_t i = 0; i < dim; ++i) {
    cout << "Computing tensor " << (double) i*(double) 100.0/(double) dim << "\% done." << endl;
    for(size_t j = 0; j < dim ; ++j) {
      for(size_t k = 0; k < dim ; ++k) {
        for(size_t l = 0; l < dim ; ++l) {
          axic[l].value().value().value()          = xic[l];
          axic[l].gradient().value().value()       = 0.0;
          axic[l].value().gradient().gradient()    = 0.0;
          axic[l].gradient().gradient().gradient() = 0.0;
          axic[l].value().value().gradient()       = 0.0;
          axic[l].gradient().value().gradient()    = 0.0;
          axic[l].value().gradient().value()       = 0.0;
          axic[l].gradient().gradient().value()    = 0.0;
        }
        axic[i].value().value().gradient() = 1.0;
        axic[j].value().gradient().value() = 1.0;
        axic[k].gradient().value().value() = 1.0;
        integrate<t3s>(axic, dim, h);
        for(size_t l = 0; l < dim ; ++l) {
          T[i][j][k][l] = axic[l].gradient().gradient().gradient();
        }
      }
      for(size_t k = 0; k < dim ; ++k) {
        H[i][j][k] = axic[k].value().gradient().gradient();
      }
    }
    for(size_t j = 0; j < dim ; ++j) {
      J[i][j] = axic[j].value().value().gradient();
    }
  }
  delete [] axic;
}

void jactest(double* xold, size_t dim, double h) {
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
      residual_beuler<t1s>(x, axold, h, y);
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
  
  
  jac_beuler<t1s>(x, axold, h, Jhc);
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
#endif
