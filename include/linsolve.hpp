#ifndef LINSOLVE_HPP
#define LINSOLVE_HPP
#include <iostream>

using namespace std;

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
      
using namespace std;

#endif
