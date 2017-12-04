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
#include "alg.hpp"
#include "user.hpp"
//#include "linsolve.hpp"

using namespace std;

using namespace alg;

typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;


void t1s_driver(double* xic, size_t dim, double h, double** J);
void t2s_t1s_driver(double* xic, size_t dim, double h, double **J, double*** H);
void t3s_t2s_t1s_driver(double* xic, size_t dim, double h,
    double **J, double*** H, double ****T);

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


void init(size_t dim) {
  py = new double[dim];
  t3_py = new double[dim];
  t1_py = new double[dim];
  t3_t1_py = new double[dim];
  t2_py = new double[dim];
  t3_t2_py = new double[dim];
  t2_t1_py = new double[dim];
  t3_t2_t1_py = new double[dim];
  pJ = new double*[dim];
  t3_pJ = new double*[dim];
  t1_pJ = new double*[dim];
  t3_t1_pJ = new double*[dim];
  t2_pJ = new double*[dim];
  t3_t2_pJ = new double*[dim];
  t2_t1_pJ = new double*[dim];
  t3_t2_t1_pJ = new double*[dim];
  pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    pJ[i] = pJ[0] + dim * i;
  }
  t3_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t3_pJ[i] = t3_pJ[0] + dim * i;
  }
  t1_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t1_pJ[i] = t1_pJ[0] + dim * i;
  }
  t3_t1_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t3_t1_pJ[i] = t3_t1_pJ[0] + dim * i;
  }
  t2_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t2_pJ[i] = t2_pJ[0] + dim * i;
  }
  t3_t2_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t3_t2_pJ[i] = t3_t2_pJ[0] + dim * i;
  }
  t2_t1_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t2_t1_pJ[i] = t2_t1_pJ[0] + dim * i;
  }
  t3_t2_t1_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    t3_t2_t1_pJ[i] = t3_t2_t1_pJ[0] + dim * i;
  }
  tmp_pJ = new double*[dim];
  tmp_pJ[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    tmp_pJ[i] = tmp_pJ[0] + dim * i;
  }
}

void destroy() {
  delete [] py; delete [] t3_py, delete [] t1_py;
  delete [] t3_t1_py;
  delete [] t2_py; delete [] t2_t1_py;
  delete [] pJ[0]; delete [] t1_pJ[0];
  delete [] t3_t2_py; delete [] t3_t2_t1_py;
  delete [] t3_pJ[0]; delete [] t3_t1_pJ[0];
  delete [] t2_pJ[0]; delete [] t2_t1_pJ[0];
  delete [] pJ; delete [] t1_pJ;
  delete [] t3_t2_pJ[0]; delete [] t3_t2_t1_pJ[0];
  delete [] t3_pJ; delete [] t3_t1_pJ;
  delete [] t2_pJ; delete [] t2_t1_pJ;
  delete [] t3_t2_pJ; delete [] t3_t2_t1_pJ;
  delete [] tmp_pJ[0];
  delete [] tmp_pJ;
}

#if 0

template <class T> int adlinsolve(T **A, T *B, size_t n) {
  T t;
  cout << "adlinsolve not implemented for this type " << endl;
}

template <> int adlinsolve<double>(double **J, double *y, size_t dim) {
  return solve(J, y, dim);
}

template <> int adlinsolve<t1s>(t1s **t1s_J, t1s *t1s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for (size_t i = 0; i < dim; ++i) {
    py[i] = t1s_y[i].getValue();
    t1_py[i] = t1s_y[i].getGradient();
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t1s_J[i][j].getValue();
      tmp_pJ[i][j] = t1s_J[i][j].getValue();
      t1_pJ[i][j] = t1s_J[i][j].getGradient();
    }
  }

  // Solve 1st order system
  int ierr = solve(pJ, py, dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t1_py, dim);
  // That's it, we have the tangents t1_x in t1_py of the LS Ax=b
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for (size_t i = 0; i < dim; ++i) {
    t1s_y[i] = py[i];
    t1s_y[i].setGradient(t1_py[i]);
  }
  return 0;
}

template <> int adlinsolve<t2s>(t2s **t2s_J, t2s *t2s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for (size_t i = 0; i < dim; ++i) {
    py[i] = t2s_y[i].value().value();
    t1_py[i] = t2s_y[i].value().gradient();
    t2_py[i] = t2s_y[i].gradient().value();
    t2_t1_py[i] = t2s_y[i].gradient().gradient();
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = t2s_J[i][j].value().value();
      tmp_pJ[i][j] = t2s_J[i][j].value().value();
      t1_pJ[i][j] = t2s_J[i][j].value().gradient();
      t2_pJ[i][j] = t2s_J[i][j].gradient().value();
      t2_t1_pJ[i][j] = t2s_J[i][j].gradient().gradient();
    }
  }

  // Solve 1st order system for t1 and t2
  int ierr = solve(pJ, py, dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  decmatmul(t2_pJ, py, t2_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t1_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t2_py, dim);
  decmatmul(t2_t1_pJ, py, t2_t1_py, dim);
  decmatmul(t1_pJ, t2_py, t2_t1_py, dim);
  decmatmul(t2_pJ, t1_py, t2_t1_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t2_t1_py, dim);
  // Put x and t1_x back into the t1s type
  // cout << "New x and step" << endl;
  for (size_t i = 0; i < dim; ++i) {
    t2s_y[i] = py[i];
    t2s_y[i].gradient().gradient() = t2_t1_py[i];
    t2s_y[i].value().gradient() = t1_py[i];
    t2s_y[i].gradient().value() = t2_py[i];
  }
  return 0;
}

template <> int adlinsolve<t3s>(t3s **t3s_J, t3s *t3s_y, size_t dim) {
  // Get the values and tangents out for both J and y
  for (size_t i = 0; i < dim; ++i) {
    py[i] = t3s_y[i].value().value().value();
    t3_py[i] = t3s_y[i].gradient().value().value();
    t1_py[i] = t3s_y[i].value().value().gradient();
    t3_t1_py[i] = t3s_y[i].gradient().value().gradient();
    t2_py[i] = t3s_y[i].value().gradient().value();
    t3_t2_py[i] = t3s_y[i].gradient().gradient().value();
    t2_t1_py[i] = t3s_y[i].value().gradient().gradient();
    t3_t2_t1_py[i] = t3s_y[i].gradient().gradient().gradient();
    for (size_t j = 0; j < dim; ++j) {
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
  int ierr = solve(pJ, py, dim);
  // t1_py has the tangents of the RHS of the primal. We now do t1_b - A_1*x
  // which is the RHS of the 1st order LS and decrement A_1*x. t1_b was already
  // extracted above
  decmatmul(t1_pJ, py, t1_py, dim);
  decmatmul(t2_pJ, py, t2_py, dim);
  decmatmul(t3_pJ, py, t3_py, dim);
  // Use the saved Jacobian. The matrix is the same for the 1st order LS
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t1_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t2_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t3_py, dim);

  decmatmul(t2_t1_pJ, py, t2_t1_py, dim);
  decmatmul(t1_pJ, t2_py, t2_t1_py, dim);
  decmatmul(t2_pJ, t1_py, t2_t1_py, dim);

  decmatmul(t3_t2_pJ, py, t3_t2_py, dim);
  decmatmul(t3_pJ, t2_py, t3_t2_py, dim);
  decmatmul(t2_pJ, t3_py, t3_t2_py, dim);

  decmatmul(t3_t1_pJ, py, t3_t1_py, dim);
  decmatmul(t1_pJ, t3_py, t3_t1_py, dim);
  decmatmul(t3_pJ, t1_py, t3_t1_py, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t2_t1_py, dim);

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t3_t2_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t3_t1_py, dim);
  decmatmul(t3_t1_pJ, t2_py, t3_t2_t1_py, dim);
  decmatmul(t1_pJ, t3_t2_py, t3_t2_t1_py, dim);
  decmatmul(t3_t2_t1_pJ, py, t3_t2_t1_py, dim);
  decmatmul(t2_t1_pJ, t3_py, t3_t2_t1_py, dim);
  decmatmul(t3_t2_pJ, t1_py, t3_t2_t1_py, dim);
  decmatmul(t2_pJ, t3_t1_py, t3_t2_t1_py, dim);
  decmatmul(t3_pJ, t2_t1_py, t3_t2_t1_py, dim);
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pJ[i][j] = tmp_pJ[i][j];
    }
  }
  ierr = solve(pJ, t3_t2_t1_py, dim);
  // Put x and t1_x back into the t1s type
  for (size_t i = 0; i < dim; ++i) {
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

/*!
   \brief "Driver for accumulating Jacobian using finite difference"
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "Jacobian"
   \pre "Input system and state"
   \post "Jacobian"
*/
void fdJ_driver(double* xic, size_t dim, double h, double** J) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double pert = 1e-8;

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    integrate<double>(xpert1, dim, h);
    integrate<double>(xpert2, dim, h);
    for (size_t j = 0; j < dim; ++j) J[i][j] = (xpert1[j] - xpert2[j])/pert;
  }
  integrate<double>(xic, dim, h);
  delete [] xpert1;
  delete [] xpert2;
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
void fdH_driver(double* xic, size_t dim, double h, double*** H) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double **Jpert1 = new double*[dim];
  Jpert1[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) Jpert1[i] = Jpert1[0] + i*dim;
  double **Jpert2 = new double*[dim];
  Jpert2[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) Jpert2[i] = Jpert2[0] + i*dim;
  double pert = 1e-8;

  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t1s_driver(xpert1, dim, h, Jpert1);
    t1s_driver(xpert2, dim, h, Jpert2);
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        H[i][j][k] = (Jpert1[j][k] - Jpert2[j][k])/pert;
      }
    }
  }
  delete [] xpert1;    delete [] xpert2;
  delete [] Jpert1[0]; delete [] Jpert1;
  delete [] Jpert2[0]; delete [] Jpert2;
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
void fdT_driver(double* xic, size_t dim, double h, double**** T) {
  double *xpert1 = new double[dim];
  double *xpert2 = new double[dim];
  double **J = new double* [dim];
  J[0] = new double[dim*dim];
  for (size_t i = 0; i < dim; ++i) J[i] = J[0] + i*dim;
  double ***Hpert1 = new double**[dim];
  Hpert1[0] = new double*[dim*dim];
  Hpert1[0][0] = new double[dim*dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    Hpert1[i] = Hpert1[0] + i * dim;
    for (size_t j = 0; j < dim; ++j) {
      Hpert1[i][j] = Hpert1[0][0] + i*dim*dim + j*dim;
    }
  }
  double ***Hpert2 = new double**[dim];
  Hpert2[0] = new double*[dim*dim];
  Hpert2[0][0] = new double[dim*dim*dim];
  for (size_t i = 0; i < dim; ++i) {
    Hpert2[i] = Hpert2[0] + i * dim;
    for (size_t j = 0; j < dim; ++j) {
      Hpert2[i][j] = Hpert2[0][0] + i*dim*dim + j*dim;
    }
  }
  double pert = 1e-8;
  for (size_t i = 0; i < dim; ++i) {
    cout << "Computing tensor. "
      << (double) i*(double) 100.0/(double) dim << "\% done." << endl;
    for (size_t j = 0; j < dim; ++j) xpert1[j] = xic[j];
    for (size_t j = 0; j < dim; ++j) xpert2[j] = xic[j];
    xpert1[i] += pert/2.0;
    xpert2[i] -= pert/2.0;
    t2s_t1s_driver(xpert1, dim, h, J, Hpert1);
    t2s_t1s_driver(xpert2, dim, h, J, Hpert2);
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        for (size_t l = 0; l < dim; ++l) {
          T[i][j][k][l] = (Hpert1[j][k][l] - Hpert2[j][k][l])/pert;
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

/*!
   \brief "Driver for accumulating Jacobian using AD. Go over all 
   Cartesian basis vectors of tangent t1_xic and collect one Jacobian column after the other."
   \param xic "Initial conditions"
   \param dim "Dimension of the state"
   \param J "Jacobian"
   \pre "Input system and state"
   \post "Jacobian"
*/
void t1s_driver(double* xic, size_t dim, double h, double** J) {
  t1s* axic = new t1s[dim];
  for (size_t i = 0; i < dim; ++i) {
    for (size_t j = 0; j < dim ; ++j) {
      axic[j] = xic[j];
      axic[j].setGradient(0.0);
    }
    axic[i].setGradient(1.0);
    integrate<t1s>(axic, dim, h);
    for (size_t j = 0; j < dim; ++j) J[i][j] = axic[j].getGradient();
  }
  for (size_t i = 0; i < dim; i++) xic[i] = axic[i].getValue();
  delete [] axic;
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
void t2s_t1s_driver(double* xic, size_t dim, double h,
    double **J, double*** H) {
  for (size_t i = 0; i < dim; ++i) {
    t2s* axic = new t2s[dim];
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        axic[k].value().value()       = xic[k];
        axic[k].gradient().gradient() = 0.0;
        axic[k].value().gradient()    = 0.0;
        axic[k].gradient().value()    = 0.0;
      }
      axic[i].value().gradient() = 1.0;
      axic[j].gradient().value() = 1.0;
      integrate<t2s>(axic, dim, h);
      for (size_t k = 0; k < dim; ++k) {
        H[i][j][k] = axic[k].gradient().gradient();
      }
      // This should give you the Jacobian too
      // for(size_t k = 0; k < dim ; ++k) {
      //   J[j][k] = axic[k].gradient().value();
      // }
    }
    for (size_t j = 0; j < dim; ++j) {
      J[i][j] = axic[j].value().gradient();
    }
    delete [] axic;
  }
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
void t3s_t2s_t1s_driver(double* xic, size_t dim, double h,
    double **J, double*** H, double ****T) {
  t3s* axic = new t3s[dim];
  for (size_t i = 0; i < dim; ++i) {
    cout << "Computing tensor " << (double) i*(double) 100.0/(double) dim
      << "\% done." << endl;
    for (size_t j = 0; j < dim; ++j) {
      for (size_t k = 0; k < dim; ++k) {
        for (size_t l = 0; l < dim; ++l) {
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
        for (size_t l = 0; l < dim; ++l) {
          T[i][j][k][l] = axic[l].gradient().gradient().gradient();
        }
      }
      for (size_t k = 0; k < dim; ++k) {
        H[i][j][k] = axic[k].value().gradient().gradient();
      }
    }
    for (size_t j = 0; j < dim; ++j) {
      J[i][j] = axic[j].value().value().gradient();
    }
  }
  delete [] axic;
}

#endif
/*!
   \brief "Test the user provided implementation of the Jacobian. Outputs the
   handwritten Jacobian and the AD generated Jacobian based on residual_beuler"
   \param xold "State xold"
*/
void jactest(double* xold, size_t dim, double h) {
  pVector<t1s> x = pVector<t1s>(dim);
  pVector<t1s> axold = pVector<t1s>(dim);
  pVector<t1s> y = pVector<t1s>(dim);
  for (size_t i = 0; i < dim; ++i) {
    axold.set(i,xold[i]);
    x.set(i,xold[i]);
  }

  // Evaluate jacobian
  pMatrix<double> J(dim, dim);
  for (size_t i = 0; i < dim; ++i) {
    x.get_datap()[i].setGradient(1.0);
    for (size_t j = 0; j < dim; ++j) {
      residual_beuler<t1s>(x, axold, h, y);
      J.set(i, j, y.get_datap()[j].getGradient());
    }
    x.get_datap()[i].setGradient(0.0);
  }

  cout << "AD Jacobian" << endl;
  // Print jacobian
  J.display();
  
  // Hand coded jacobian

  pMatrix<double> Jhc(dim, dim);
  pVector<double> xold_hc(dim);
  pVector<double> x_hc(dim);
  for (size_t i = 0; i < dim; ++i) {
    xold_hc.set(i,xold[i]);
    x_hc.set(i,xold[i]);
  }

  for (size_t i = 0; i < dim; ++i) xold_hc.set(i, xold[i]);

  cout << "HC Jacobian" << endl;
  jac_beuler<double>(x_hc, xold_hc, h, Jhc);

  Jhc.display();

}
#endif  // ADUPROP_AD_HPP_
