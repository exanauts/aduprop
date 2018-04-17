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
#include <fstream>
// AD library used
// Rudimentary linear solver interface
#include "ad.hpp"
//#include "cxxopts.hpp"
#include "parallel.hpp"

int main(int argc, char* argv[]) {

  // activate timer in propagateAD
  
  global_prof.activate("propagateAD");
  global_prof.activate("propagateCOV");
  global_prof.activate("propagateH");
  global_prof.activate("propagateMU");
  global_prof.activate("reduction");

  // Variable declaration
  size_t tsteps = 1;
  // pVector<double> xold, x, F;
  // pVector<double> x0;
  // pMatrix<double> TMAT;
  System sys;
  size_t dim = sys.dim();

  pVector<double> xold(dim);
  pVector<double> x(dim);
  // F.alloc(sys.dimension);

  // Read initial conditions from external file. Check for dimension
  // consistency.
  
  ad drivers(sys);

  dim = sys.dimension;

  // Co-variance
  pMatrix<double> cv0(sys.dimension, sys.dimension);
  cv0.zeros();
  // Strings
  std::string fcov;

  // Where do we put this???
  // We should only alocate if we're using all of these.
  size_t chunk = paduprop_getend(dim) - paduprop_getstart(dim);
  //pTensor4<double> T(dim, dim, dim, chunk);
  //pTensor3<double> H(dim, dim, dim);
  pTensor4<double> T;
  pTensor3<double> H(dim, dim, chunk);
  pMatrix<double>  J(dim, dim);

  for (size_t i = 0; i < sys.dimension; ++i) 
    cv0[i][i] = 0.0000001;
  
  cv0[sys.dim()][sys.dimension] = 0.000001;

  for (size_t i = 0; i < tsteps; ++i) {
    // Save mean in trajectory matrix
    // for (size_t j = 0; j < sys.dimension; ++j) {
    //   TMAT.set(j, i, x[j]);
    // }

    // Propagate
    std::cout << "Step: " << i << "." << std::endl;
    int degree = 2;
    propagateAD(x, cv0, sys, J, H, T, drivers, degree);
  }
  if(paduprop_getrank() == 0) {  
    double compare, diff;
    std::ifstream infile;
    infile.open("solution.txt");
    infile >> compare;
    diff = compare - cv0.norm();
    std::cout << std::setprecision(16) << "Sol: " << compare << ". " << std::endl;
    std::cout << std::setprecision(16) << "cv0 norm: " << cv0.norm() << "." << std::endl;
    std::cout << std::setprecision(16) << "Diff: " << diff << "." << std::endl;
    infile.close();
  }

  std::cout << global_prof << std::endl;

  return 0;
}

