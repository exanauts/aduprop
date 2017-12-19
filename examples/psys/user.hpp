/*!
   \file "user.hpp"
   \brief "This is the provided user code."
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/
#ifndef USER_HPP
#define USER_HPP
#include <iostream>
#include "alg.hpp"

typedef class System {
public:
  // System parameters.
  // Generator
  
  double h; 
  size_t dimension;
  

size_t dim() { return dimension; }

/*!
   \brief "Generic residual function written by the user. All variables that
   are on the computational path from an independent to dependent need to be of
   type T. All other variables may be passive (double)."
   \param x "New x"
   \param xold "Old x"
   \param F "Residual"
   \param h "Discretization step"
   \pre "System new state x and old state xold"
   \post "Residual F"
*/
template <class T> void residual_beuler(const alg::pVector<T> &x,
    const alg::pVector<T> &xold, alg::pVector<T> &F) {
  // Empty
}

/*!
   \brief "User provided Jacobian"
   \param x "New x"
   \param xold "Old x"
   \param J "Jacobian"
   \param h "Discretization step"
   \pre "Input states x and xold"
   \post "Jacobian J"
*/
template <class T> void jac_beuler(const alg::pVector<T> &x, 
    const alg::pVector<T> &xold, alg::pMatrix<T> &J) {
  // Empty
}

} System;
#endif
