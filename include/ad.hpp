/*!
\file "ad.hpp"

\author "Adrian Maldonado and Michel Schanen"

\date 21/11/2017

\brief Driver for Jacobian, Hessian and 3rd order tensor provided by
Automatic Differentiation. The drivers call the user provided residual, Jacobian
and integration.

We use the AD notation from "The Art of Differentiating Computer
Programs", by Uwe Naumann.

y = f(x)
1st order tangent-linear:
y = f(x)
y_1 = <f'(x), x_1> 
2nd order tangent-linear:
y = f(x)
y_2 = <f'(x), x_2> 
y_1 = <f'(x), x_1> 
y_(1,2) = <f''(x), x_1, x_2> + <f'(x), x_(1,2)> 
3rd order tangent-linear has 8 statements. There will be 8 variables for every
original variable.
*/

#ifndef ADUPROP_AD_HPP_
#define ADUPROP_AD_HPP_
#include <iostream>
#include <cstdlib> 
#include "perf.hpp"
#include "alg.hpp"
#include "user.hpp"
#include <stdlib.h>
#include "parallel.hpp"
#ifdef EIGEN
#include <Eigen/Dense>
#endif
#ifdef EIGEN_SPARSE
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU> 
#endif
#ifdef __INTEL_COMPILER
#define RESTRICT restrict
#else
#define RESTRICT __restrict__
#endif
    

using namespace std;
using namespace alg;
#ifdef EIGEN
using namespace Eigen;
#endif
#ifdef EIGEN_SPARSE
using namespace Eigen;
#endif

typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

perf global_prof("global");

extern System sys;

// TODO: Name is not informative. This class does much more than AD: stores data, has integration functions,
// etc...
typedef class ad {
private:
  pVector<double> py;
  pVector<double> t3_py;
  pVector<double> t1_py;
  pVector<double> t3_t1_py;
  pVector<double> t2_py;
  pVector<double> t3_t2_py;
  pVector<double> t2_t1_py;
  pVector<double> t3_t2_t1_py;
  pMatrix<double> pJ;
  pMatrix<double> t3_pJ;
  pMatrix<double> t1_pJ;
  pMatrix<double> t3_t1_pJ;
  pMatrix<double> t2_pJ;
  pMatrix<double> t3_t2_pJ;
  pMatrix<double> t2_t1_pJ;
  pMatrix<double> t3_t2_t1_pJ;
  pMatrix<double> tmp_pJ;
  System *sys;
  perf prof;

public:
ad(System &sys_) : prof("AD") {
  sys = &sys_;
  size_t dim = sys->dim();
  py = pVector<double>(dim);
  t3_py = pVector<double>(dim);
  t1_py = pVector<double>(dim);
  t3_t1_py = pVector<double>(dim);
  t2_py = pVector<double>(dim);
  t3_t2_py = pVector<double>(dim);
  t2_t1_py = pVector<double>(dim);
  t3_t2_t1_py = pVector<double>(dim);
  pJ = pMatrix<double>(dim, dim);
  t3_pJ = pMatrix<double>(dim, dim);
  t1_pJ = pMatrix<double>(dim, dim);
  t3_t1_pJ = pMatrix<double>(dim, dim);
  t2_pJ = pMatrix<double>(dim, dim);
  t3_t2_pJ = pMatrix<double>(dim, dim);
  t2_t1_pJ = pMatrix<double>(dim, dim);
  t3_t2_t1_pJ = pMatrix<double>(dim, dim);
  tmp_pJ = pMatrix<double>(dim, dim);
  prof.activate("integrate");
  prof.activate("driver");
  prof.activate("t1s_driver");
  prof.activate("t2s_t1s_driver");
  prof.activate("t3s_t2s_t1s_driver");
  prof.activate("adlinsolve");
  prof.activate("t1s_adlinsolve");
  prof.activate("t2s_adlinsolve");
  prof.activate("t3s_adlinsolve");
  paduprop_init();
  
}

~ad() {
  std::cout << prof << std::endl;
  paduprop_destroy();
}

template <class T> void adlinsolve(pMatrix<T> &J, pVector<T> &y) {
  T t;
  cout << "adlinsolve not implemented for this type " << endl;
}


/*!
   \brief "Driver for accumulating Jacobian using finite difference"
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "Jacobian"
   \pre "Input system and state"
   \post "Jacobian"
*/
void fdJ_driver(const pVector<double> &xic, pMatrix<double> &J) {
  size_t dim = xic.dim();
  pVector<double> xpert1(dim);
  pVector<double> xpert2(dim);
  double pert = 1e-10;

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    integrate<double>(xpert1);
    integrate<double>(xpert2);
    for (size_t j = 0; j < dim; ++j) J[j][i] = (xpert1[j] - xpert2[j])/pert;
  }
}

/*!
   \brief "Driver for accumulating Hessian using finite difference. To avoid
   numerical issues it relies on Jacobian generated via AD"
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param H "Hessian"
   \pre "Input system and state"
   \post "Hessian"
*/
void fdH_driver(const pVector<double> &xic, pTensor3<double> &H, size_t start = 0, size_t end =0) {
  size_t dim = xic.dim();
  // no end argument is set, hence pick dim as end
  if(end == 0) end = dim;
  pVector<double> xpert1(dim);
  pVector<double> xpert2(dim);
  pMatrix<double> Jpert1(dim, dim);
  pMatrix<double> Jpert2(dim, dim);
  double pert = 1e-8;
  for (size_t i = start; i < end; ++i) {
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t1s_driver(xpert1, Jpert1);
    t1s_driver(xpert2, Jpert2);
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        H[j][k][i-start] = (Jpert1[j][k] - Jpert2[j][k])/pert;
      }
    }
  }
}

/*!
   \brief "Driver for accumulating 3rd order tensor using finite difference. To avoid
   numerical issues it relies on Hessian generated via AD"
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param sys "System parameters"
   \param T "Tensor"
   \pre "Input system and state"
   \post "Tensor"
*/
void fdT_driver(const pVector<double> &xic, pTensor4<double> &T, size_t start = 0, size_t end = 0) {
  size_t dim = xic.dim();
  if(end == 0) end = dim;
  pVector<double> xpert1(dim);
  pVector<double> xpert2(dim);
  pMatrix<double> J(dim, dim);
  pTensor3<double> Hpert1(dim, dim, dim);
  pTensor3<double> Hpert2(dim, dim, dim);
  double pert = 1e-8;
  for (size_t i = start; i < end; ++i) {
    cout << "Computing tensor. "
      << (double) i*(double) 100.0/(double) dim << "\% done." << endl;
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t2s_t1s_driver(xpert1, J, Hpert1);
    t2s_t1s_driver(xpert2, J, Hpert2);
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        for (size_t l = 0; l < dim; ++l) {
          T[j][k][l][i-start] = (Hpert1[j][k][l] - Hpert2[j][k][l])/pert;
        }
      }
    }
  }
}

/*!
   \brief "Driver for accumulating Jacobian using AD. Go over all 
   Cartesian basis vectors of tangent t1_xic and collect one Jacobian column after the other."
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "Jacobian"
   \pre "Input system and state"
   \post "Jacobian"
*/
void t1s_driver(const pVector<double> &xic, pMatrix<double> &J) {
  prof.begin("t1s_driver");
  size_t dim = xic.dim();
  pVector<t1s> axic(dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim ; ++j) {
      axic[j] = xic[j];
      axic[j].setGradient(0.0);
    }
    axic[i].setGradient(1.0);
    integrate<t1s>(axic);
    for (size_t j = 0; j < dim; ++j) J[j][i] = axic[j].getGradient();
  }
  prof.end("t1s_driver");
}


/*!
   \brief "Driver for accumulating Hessian using AD. Go over all 
   Cartesian basis vectors of the tangents t1_xic and t2_xic and collect one Hessian projection after the other."
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "Hessian"
   \pre "Input system and state"
   \post "Hessian"
*/
void t2s_t1s_driver(const pVector<double> &xic,
    pMatrix<double> &J, pTensor3<double> &H, size_t start = 0, size_t end = 0) {
  prof.begin("t2s_t1s_driver");
  size_t dim = xic.dim();
  
  // no end argument is set, hence pick dim as end
  if(end == 0) end = dim;
  pVector<t2s> axic(dim);
  for (size_t i = start; i < end; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        axic[k].value().value()       = xic[k];
        axic[k].gradient().gradient() = 0.0;
        axic[k].value().gradient()    = 0.0;
        axic[k].gradient().value()    = 0.0;
      }
      axic[i].value().gradient() = 1.0;
      axic[j].gradient().value() = 1.0;
      integrate<t2s>(axic);
      for (size_t k = 0; k < dim; ++k) {
        H[k][j][i-start] = axic[k].gradient().gradient();
      }
    }
  }
  prof.end("t2s_t1s_driver");
}

/*!
   \brief "Driver for accumulating 3rd order tensor using AD. Go over all 
   Cartesian basis vectors of the tangents t1_xic, t2_xic and t3_xic 
   and collect one tensor projection after the other."
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "tensor"
   \pre "Input system and state"
   \post "tensor"
*/
void t3s_t2s_t1s_driver(const pVector<double> &xic,
    pMatrix<double> &J, pTensor3<double> &H, pTensor4<double> &T, size_t start = 0, size_t end = 0) {
  prof.begin("t3s_t2s_t1s_driver");
  size_t dim = xic.dim();
  
  // no end argument is set, hence pick dim as end
  if(end == 0) end = dim;
  pVector<t3s> axic(dim);
  for (size_t i = start; i < end; ++i) {
    if(paduprop_getrank() == 0) {
    cout << "Computing tensor " << (double) (i-start)*(double) 100.0/(double) (end-start)
      << "\% done." << endl;
    }
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k)
      {
        for (size_t l = 0; l < dim; ++l)
        {
          axic[l].value().value().value() = xic[l];
          axic[l].gradient().value().value() = 0.0;
          axic[l].value().gradient().gradient() = 0.0;
          axic[l].gradient().gradient().gradient() = 0.0;
          axic[l].value().value().gradient() = 0.0;
          axic[l].gradient().value().gradient() = 0.0;
          axic[l].value().gradient().value() = 0.0;
          axic[l].gradient().gradient().value() = 0.0;
        }
        axic[i].value().value().gradient() = 1.0;
        axic[j].value().gradient().value() = 1.0;
        axic[k].gradient().value().value() = 1.0;
        integrate<t3s>(axic);
        for (size_t l = 0; l < dim; ++l)
        {
          T[l][k][j][i - start] = axic[l].gradient().gradient().gradient();
        }
      }
    }
  }
  prof.end("t3s_t2s_t1s_driver");
}

/*!
   \brief "Test the user provided implementation of the Jacobian. Outputs the
   handwritten Jacobian and the AD generated Jacobian based on residual_beuler"
   \param xold "State xold"
*/
void jactest(const pVector<double> &xold) {
  size_t dim = xold.dim();
  pVector<t1s> x = pVector<t1s>(dim);
  pVector<t1s> axold = pVector<t1s>(dim);
  pVector<t1s> y = pVector<t1s>(dim);
  for (size_t i = 0; i < dim; ++i) {
    axold[i] = xold[i];
    x[i] = xold[i];
  }

  // Evaluate jacobian
  pMatrix<double> J(dim, dim);
  for (size_t j = 0; j < dim; ++j) {
    x[j].setGradient(1.0);
    for (size_t i = 0; i < dim; ++i) {
      sys->residual_beuler<t1s>(x, axold, y);
      J[i][j] = y[i].getGradient();
    }
    x[j].setGradient(0.0);
  }

  cout << "AD Jacobian" << endl;
  // Print jacobian
  cout << J;

  // Hand coded jacobian

  pMatrix<double> Jhc(dim, dim);
  pVector<double> xold_hc(dim);
  pVector<double> x_hc(dim);
  for (size_t i = 0; i < dim; ++i) {
    xold_hc[i] = xold[i];
    x_hc[i] = xold[i];
    for (size_t j = 0; j < dim; ++j) {
      Jhc[i][j] = 0;
    }
  }

  cout << "HC Jacobian" << endl;
  sys->jac_beuler<double>(x_hc, xold_hc, Jhc);
  cout << Jhc;
  double diff = (J - Jhc).norm();
  cout << "Norm" << endl << diff << endl;
  if(diff > 1e-14) {
    cerr << "ERROR: HC and AD Jacobian differ." << endl;
    exit(1);
  }
  
}

/*!
   \brief "User provided time integration"
   \param x "1st order active input variable"
   \param h "Discretization"
   \pre "Initial conditions with input tangents"
   \post "New state x with 1st order tangents"
*/
template <class T> void integrate(pVector<T> &x) { 
  // TODO(Michel, Adrian): integrate must not overwrite &x by default.

  prof.begin("integrate");
  
  size_t dim = x.dim();
  double eps = 1e-9;
  int iteration = 0;
  pVector<T> xold(dim);
  pVector<T> yold(dim);
  pVector<T> y(dim);
  pMatrix<T> J(dim, dim);
  pMatrix<T> Jold(dim, dim);
  xold = x;
  sys->residual_beuler<T>(x, xold, y);
  J.zeros();
#ifdef EIGEN
  Map<Matrix<T, Dynamic, 1> > eigymap(y.get_datap(), y.dim());
  Map<Matrix<T, Dynamic, Dynamic> > eigJmap(J.get_datap(), J.nrows(), J.ncols());
#endif
#ifdef EIGEN_SPARSE
  Map<Matrix<T, Dynamic, 1> > eigymap(y.get_datap(), y.dim());
  SparseMatrix<T> sJ(y.dim(), y.dim());
  std::vector<Eigen::Triplet<T> > tripletList;
  tripletList.reserve(y.dim()*y.dim());
  SparseLU<SparseMatrix<T> >  solver;
#ifdef EIGEN_SPARSE
  sys->jac_beuler<T>(x, xold, J);
  for(size_t i = 0; i < J.nrows(); ++i) {
    for (size_t j = 0; j < J.ncols(); ++j) {
      if(J[i][j] != 0.0) {
        tripletList.push_back(Eigen::Triplet<T>(i,j, J[i][j]));
      }
    }
  }
  sJ.setFromTriplets(tripletList.begin(), tripletList.end());
  tripletList.clear();
  solver.analyzePattern(sJ); 
#endif
#endif
  do {
    iteration = iteration + 1;
    sys->jac_beuler<T>(x, xold, J);
#ifdef DBUG
    cout << endl << endl << "Iteration: " << iteration << endl;
    cout << "Norm(J): " << setprecision(20) << J.norm() << endl << " J: " << J << endl;
    cout << "Norm(F): " << y.norm() << endl << " F: " << y << endl;
#endif
    yold = y;
    Jold = J;
#ifdef EIGEN_SPARSE
    for(size_t i = 0; i < J.nrows(); ++i) {
      for (size_t j = 0; j < J.ncols(); ++j) {
        if(J[i][j] != 0.0) {
          tripletList.push_back(Eigen::Triplet<T>(i,j, J[i][j]));
        }
      }
    }
    sJ.setFromTriplets(tripletList.begin(), tripletList.end());
    tripletList.clear();
    solver.compute(sJ);
    eigymap = solver.solve(eigymap);
#elif EIGEN
    eigymap = eigJmap.fullPivLu().solve(eigymap);
#else
    adlinsolve<T>(J, y);
#endif
    pVector<T> res(dim);
    res = Jold * y - yold;
#ifdef DBUG
    cout << "Norm(y): " << y.norm() << endl << " y: " << y << endl;
    cout << "|Ax - b| " << res.norm() << endl;
    cout << "|Ax - b| " << res << endl;
#endif
    x = x - y;
    sys->residual_beuler<T>(x, xold, y);
  } while (y.norm() > eps);
  
  prof.end("integrate");
}

} ad;

template <> void ad::adlinsolve<double>(pMatrix<double> &J, 
    pVector<double> &y) {
  prof.begin("adlinsolve");
  LUsolve(J, y);
  prof.end("adlinsolve");
}


template <> void ad::adlinsolve<t1s>(pMatrix<t1s> &t1s_J, pVector<t1s> &t1s_y) {
  // Get the values and tangents out for both J and y
  prof.begin("t1s_adlinsolve");
  size_t dim = t1s_y.dim();
  for (size_t i = 0; i < dim; ++i) py[i] = t1s_y[i].getValue();
  for (size_t i = 0; i < dim; ++i) t1_py[i] = t1s_y[i].getGradient();
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t1s_J[i][j].getValue();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      tmp_pJ[i][j] = t1s_J[i][j].getValue();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t1_pJ[i][j] = t1s_J[i][j].getGradient();
    }
  }

  // Solve 1st order system
  LUsolve(pJ, py);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  decmatmul(t1_pJ, py, t1_py);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t1_py);
  // That's it, we have the tangents t1_x in t1_py of the LS Ax=b
  // Put x and t1_x back into the t1s type
  for (size_t i = 0; i < dim; ++i) {
    t1s_y[i] = py[i];
    t1s_y[i].setGradient(t1_py[i]);
  }
  prof.end("t1s_adlinsolve");
}


template <> void ad::adlinsolve<t2s>(pMatrix<t2s> &t2s_J, pVector<t2s> &t2s_y) {
  // Get the values and tangents out for both J and y
  prof.begin("t2s_adlinsolve");
  size_t dim = t2s_y.dim();

  for (size_t i = 0; i < dim; ++i) py[i] = t2s_y[i].value().value();
  for (size_t i = 0; i < dim; ++i) t1_py[i] = t2s_y[i].value().gradient();
  for (size_t i = 0; i < dim; ++i) t2_py[i] = t2s_y[i].gradient().value();
  for (size_t i = 0; i < dim; ++i) t2_t1_py[i] = t2s_y[i].gradient().gradient();
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t2s_J[i][j].value().value();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      tmp_pJ[i][j] = t2s_J[i][j].value().value();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t1_pJ[i][j] = t2s_J[i][j].value().gradient();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t2_pJ[i][j] = t2s_J[i][j].gradient().value();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t2_t1_pJ[i][j] = t2s_J[i][j].gradient().gradient();
    }
  }

  // Solve 1st order system for t1 and t2
  LUsolve(pJ, py);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  decmatmul(t1_pJ, py, t1_py);
  decmatmul(t2_pJ, py, t2_py);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t2_py);
  decmatmul(t2_t1_pJ, py, t2_t1_py);
  decmatmul(t1_pJ, t2_py, t2_t1_py);
  decmatmul(t2_pJ, t1_py, t2_t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t2_t1_py);
  // Put x and t1_x back into the t1s type

  for (size_t i = 0; i < dim; ++i) t2s_y[i] = py[i];
  for (size_t i = 0; i < dim; ++i) t2s_y[i].gradient().gradient() = t2_t1_py[i];
  for (size_t i = 0; i < dim; ++i) t2s_y[i].value().gradient() = t1_py[i];
  for (size_t i = 0; i < dim; ++i) t2s_y[i].gradient().value() = t2_py[i];
  prof.end("t2s_adlinsolve");
}

template <> void ad::adlinsolve<t3s>(pMatrix<t3s> &t3s_J, pVector<t3s> &t3s_y) {
  prof.begin("t3s_adlinsolve");
  // Get the values and tangents out for both J and y
  size_t dim = t3s_y.dim();
  for (size_t i = 0; i < dim; ++i) py[i] = t3s_y[i].value().value().value();
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t3s_J[i][j].value().value().value();
    }
  }
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      tmp_pJ[i][j] = t3s_J[i][j].value().value().value();
    }
  }

  // Solve 1st order system for t1 and t2
  LUsolve(pJ, py);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t1_pJ[i][j] = t3s_J[i][j].value().value().gradient();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t1_py[i] = t3s_y[i].value().value().gradient();

  decmatmul(t1_pJ, py, t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t2_pJ[i][j] = t3s_J[i][j].value().gradient().value();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t2_py[i] = t3s_y[i].value().gradient().value();

  decmatmul(t2_pJ, py, t2_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t3_pJ[i][j] = t3s_J[i][j].gradient().value().value();
    }
  }

  for (size_t i = 0; i < dim; ++i)
    t3_py[i] = t3s_y[i].gradient().value().value();

  decmatmul(t3_pJ, py, t3_py);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t2_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t3_py);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t2_t1_pJ[i][j] = t3s_J[i][j].value().gradient().gradient();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t2_t1_py[i] = t3s_y[i].value().gradient().gradient();

  decmatmul(t2_t1_pJ, py, t2_t1_py);
  decmatmul(t1_pJ, t2_py, t2_t1_py);
  decmatmul(t2_pJ, t1_py, t2_t1_py);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t3_t2_pJ[i][j] = t3s_J[i][j].gradient().gradient().value();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t3_t2_py[i] = t3s_y[i].gradient().gradient().value();

  decmatmul(t3_t2_pJ, py, t3_t2_py);
  decmatmul(t3_pJ, t2_py, t3_t2_py);
  decmatmul(t2_pJ, t3_py, t3_t2_py);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t3_t1_pJ[i][j] = t3s_J[i][j].gradient().value().gradient();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t3_t1_py[i] = t3s_y[i].gradient().value().gradient();

  decmatmul(t3_t1_pJ, py, t3_t1_py);
  decmatmul(t1_pJ, t3_py, t3_t1_py);
  decmatmul(t3_pJ, t1_py, t3_t1_py);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t2_t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t3_t2_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t3_t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      t3_t2_t1_pJ[i][j] = t3s_J[i][j].gradient().gradient().gradient();
    }
  }
  for (size_t i = 0; i < dim; ++i)
    t3_t2_t1_py[i] = t3s_y[i].gradient().gradient().gradient();

  decmatmul(t3_t1_pJ, t2_py, t3_t2_t1_py);
  decmatmul(t1_pJ, t3_t2_py, t3_t2_t1_py);
  decmatmul(t3_t2_t1_pJ, py, t3_t2_t1_py);
  decmatmul(t2_t1_pJ, t3_py, t3_t2_t1_py);
  decmatmul(t3_t2_pJ, t1_py, t3_t2_t1_py);
  decmatmul(t2_pJ, t3_t1_py, t3_t2_t1_py);
  decmatmul(t3_pJ, t2_t1_py, t3_t2_t1_py);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  LUsolve(pJ, t3_t2_t1_py);
  // Put x and t1_x back into the t1s type
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].value().value().value() = py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].gradient().value().value() = t3_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].value().gradient().gradient() = t2_t1_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].gradient().gradient().gradient() = t3_t2_t1_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].value().value().gradient() = t1_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].gradient().value().gradient() = t3_t1_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].value().gradient().value() = t2_py[i];
  for (size_t i = 0; i < dim; ++i)
    t3s_y[i].gradient().gradient().value() = t3_t2_py[i];
  prof.end("t3s_adlinsolve");
}


/*!
   \brief "Propagates the mean and the covariance matrices one step."
   \param m0 "Mean"
   \param cv0 "Covariance matrix"
   \param sys "System data"
   \param drivers "AD drivers"
*/
void propagateAD(pVector<double>& m0, pMatrix<double>& cv0, System& sys,
    pMatrix<double>& J, pTensor3<double>& H, pTensor4<double>& T,
    ad& drivers, int degree = 1) {
  global_prof.begin("propagateAD");
  size_t dim = sys.dim();
  size_t start = paduprop_getstart(dim);
  size_t end = paduprop_getend(dim);
  size_t chunk = end-start;
 
  // THIS WILL BE AN ENUM THAT WE PASS
  // 1: Jacobian
  // 2: Jacobian + Hessian
  // 3: Jacobian + Hessian + Tensor-3
  // In this case we do not allocate for all elements.
  
  pVector<double>  m0_temp(dim);
  m0_temp.zeros();

  // CV_TEMP
  pMatrix<double>  cv_temp(dim, dim);
  cv_temp.zeros();
  
  pMatrix<double>  cv_temp2(dim, dim);
  cv_temp2.resize(dim, dim);
  cv_temp2.zeros();
  
  pTensor3<double> H_tmp;
  
  switch(degree) {
    case 3:
      T.zeros();
    case 2:
      H.zeros();
    case 1:
      J.zeros();
      break;
    default:
      std::cout << "Invalid option" << std::endl;
      exit(-1);
  }

  // Obtain tensors
  if(paduprop_getrank() == 0) {
    std::cout << "Obtaining tensors" << std::endl;
  }
  global_prof.begin("propagateH");
  switch(degree) {
    case 3:
      drivers.t2s_t1s_driver(m0, J, H);
      drivers.t3s_t2s_t1s_driver(m0, J, H, T, start, end);
      drivers.t1s_driver(m0, J);
      break;
    case 2:
      drivers.t2s_t1s_driver(m0, J, H, start, end);
      drivers.t1s_driver(m0, J);
      break;
    case 1:
      drivers.t1s_driver(m0, J);
      break;
    default:
      std::cout << "Invalid option" << std::endl;
      exit(-1);
  }
  global_prof.end("propagateH");
  if(paduprop_getrank() == 0) {
    std::cout << "Done with tensors" << std::endl;
  }
  
  if(paduprop_getrank() == 0) {
    std::cout << "Propagating mean" << std::endl;
  }
  // Propagate mean
  global_prof.begin("propagateMU");
  drivers.integrate(m0);

  if (degree == 2) {
    for (size_t p = 0; p < dim; ++p) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = start; j < end; ++j) {
          m0_temp[p] += (1.0/2.0)*H[p][i][j-start]*cv0[i][j];
        }
      }
    }
    paduprop_sum(m0_temp);
    m0 = m0 + m0_temp;
  }
  
  if (degree > 2) {
    for (size_t p = 0; p < dim; ++p) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
          m0[p] += (1.0/2.0)*H[p][i][j]*cv0[i][j];
        }
      }
    }
  }
  global_prof.end("propagateMU");

  if(paduprop_getrank() == 0) {
    std::cout << "Propagating covariance" << std::endl;
  }

  // Propagate covariance
  global_prof.begin("propagateCOV");

  for (size_t pm = start; pm < end; ++pm) {
    for (size_t pn = 0; pn < dim; ++pn) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
          cv_temp[pn][pm] += 0.5*((J[pn][j]*J[pm][i] +
                J[pm][j]*J[pn][i])*cv0[i][j]);
        }
      }
    }
  }
  paduprop_gather(cv_temp);

  if (degree > 2) {

    double aux, kurt;

    double * RESTRICT ptr_cv0 = cv0.get_datap();
    double * RESTRICT ptr_cv_temp2 = cv_temp2.get_datap();
    double * RESTRICT ptr_J = J.get_datap();
    double * RESTRICT ptr_H = H.get_datap();
    double * RESTRICT ptr_T = T.get_datap();

    for (size_t pm = 0; pm < dim; ++pm) { 
      for (size_t pn = 0; pn < dim; ++pn) { 
        for (size_t i = 0; i < dim; ++i) { 
          for (size_t j = 0; j < dim; ++j) { 
            for (size_t l = start; l < end; ++l) { 
              for (size_t k = 0; k < dim; ++k) { 
                aux = 0;
                kurt = ptr_cv0[i+j*dim]*ptr_cv0[k+l*dim] + ptr_cv0[i+l*dim]*ptr_cv0[j+k*dim] + ptr_cv0[i+k*dim]*ptr_cv0[l+j*dim];
                if(kurt != 0) {
                  aux += ptr_J[pn+i*dim]*ptr_T[pm*dim*dim*chunk+j*dim*chunk+k*chunk+l-start];
                  aux += ptr_J[pn+j*dim]*ptr_T[pm*dim*dim*chunk+i*dim*chunk+k*chunk+l-start];
                  aux += ptr_J[pn+k*dim]*ptr_T[pm*dim*dim*chunk+i*dim*chunk+j*chunk+l-start];
                  aux += ptr_H[pn*dim*dim+i*dim+j]*ptr_H[pm*dim*dim+k*dim+l];
                  aux += ptr_H[pn*dim*dim+i*dim+k]*ptr_H[pm*dim*dim+j*dim+l];
                  aux += ptr_H[pn*dim*dim+i*dim+l]*ptr_H[pm*dim*dim+j*dim+k];
                  aux += ptr_H[pn*dim*dim+j*dim+k]*ptr_H[pm*dim*dim+i*dim+l];
                  aux += ptr_H[pn*dim*dim+j*dim+l]*ptr_H[pm*dim*dim+i*dim+k];
                  aux += ptr_H[pn*dim*dim+k*dim+l]*ptr_H[pm*dim*dim+i*dim+j];
                  aux += ptr_T[pn*dim*dim*chunk+j*dim*chunk+k*chunk+l-start]*ptr_J[pm+i*dim];
                  aux += ptr_T[pn*dim*dim*chunk+i*dim*chunk+k*chunk+l-start]*ptr_J[pm+j*dim];
                  aux += ptr_T[pn*dim*dim*chunk+i*dim*chunk+j*chunk+l-start]*ptr_J[pm+k*dim];

                  aux *= (1.0/(24.0))*kurt;
                }
                ptr_cv_temp2[pn+dim*pm] += aux;
              }
            }
          }
        }
      }
    }
    if(paduprop_getrank() == 0) {
      std::cout << "Propagating covariance 3rd order part2" << std::endl;
    }
    for (size_t pn = 0; pn < dim; ++pn) { 
      for (size_t pm = 0; pm < dim; ++pm) { 
        for (size_t i = 0; i < dim; ++i) { 
          for (size_t j = 0; j < dim; ++j) { 
            for (size_t k = start; k < end; ++k) { 
              for (size_t l = 0; l < dim; ++l) { 
                aux = 0;
                kurt = ptr_cv0[i+dim*j]*ptr_cv0[k+dim*l] + ptr_cv0[i+dim*l]*ptr_cv0[j+dim*k] + ptr_cv0[i+dim*k]*ptr_cv0[l+dim*j];
                if(kurt != 0) {
                  aux += ptr_T[pn*dim*dim*chunk+i*dim*chunk+j*chunk+k-start]*ptr_J[pm+dim*l];
                  aux += ptr_J[pn+dim*l]*ptr_T[pm*dim*dim*chunk+i*dim*chunk+j*chunk+k-start];
                  aux *= (1.0/(24.0))*kurt;
                  aux -= (1.0/2.0)*(ptr_H[pn*dim*dim+i*dim+j]*ptr_H[pm*dim*dim+k*dim+l])*ptr_cv0[i+dim*j]*ptr_cv0[k+dim*l];
                }
                ptr_cv_temp2[pn+dim*pm] += aux;
              }
            }
          }
        }
      }
    }
    if(paduprop_getrank() == 0) {
      std::cout << "Done with covariance" << std::endl;
    }
    global_prof.begin("reduction");
    if(paduprop_getrank() == 0) {
      std::cout << "Start reduction" << std::endl;
    }
    paduprop_sum(cv_temp2);
    if(paduprop_getrank() == 0) {
      std::cout << "Done with reduction" << std::endl;
    }

    global_prof.end("reduction");
  }

  cv0 = cv_temp + cv_temp2;
  // Here would be the cutoff of the covariance matrix according to a 
  // 0 <= cutoff rate <= 1 set before.
  // Example: cv0.cutoff(cutrate);
  global_prof.end("propagateCOV");
  global_prof.end("propagateAD");
}



#endif  // ADUPROP_AD_HPP_
