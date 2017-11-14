#ifndef LINSOLVE_HPP
#define LINSOLVE_HPP

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
      
// Decremental matmul
void decmatmul(double **A, double *x, double *y, size_t n) {
  // for(size_t i = 0; i < n; ++i) {
  //   for(size_t j = 0; j < n; ++j) {
  //     y[i] -= A[j][i] * x[j];
  //   }
  // }
  int n_=n;
  char trans = 'N';
  double alpha = -1.0;
  double beta = 1.0;
  int incx=1;
  int incy=1;
  FNAME(dgemv)(&trans, &n_, &n_, &alpha, &A[0][0], &n_, x, &incx, &beta, y, &incy);
}
      
int solve(double **A, double *B, int n) {
  int lda=n;
  int ldb=n;
  int nrhs=1;
  int info=0;
  int ipiv[n];
  FNAME(dgesv)(&n, &nrhs, &A[0][0], &lda, ipiv, &B[0], &ldb, &info);
  return info;
}
#endif
