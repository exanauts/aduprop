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
#include "cxxopts.hpp"

int main(int argc, char* argv[]) {

  // activate timer in propagateAD
  
  global_prof.activate("propagateAD");

  // Variable declaration
  int nbuses, nbranches, ngen, nload;
  int tsteps = 20;
  size_t dim;
  pVector<double> xold, x, F;
  pVector<double> x0;
  pMatrix<double> TMAT;
  System sys;

  // problem definition

  nbuses = 2;
  nbranches = 1;
  ngen = 1;
  nload = 1;

  sys.init(nbuses, nbranches, ngen, nload);

  sys.branches[0].set(0, 1, 0.0001, 0.0576, 0.000);
  sys.loads[0].set(1, 1.0648453, 0.38835684);
  sys.gens[0].set(0, 1.575, 1.512, 0.29, 0.39, 0.1733,
      0.1733, 0.0787, 3.38, 6.09, 1.0, 0.05, 0.15);
  sys.gens[0].e_fd = 2.36980307364616349375;
  sys.gens[0].p_m = 1.06496000000000012875;

  sys.build_ybus();

  xold.alloc(sys.dimension);
  x.alloc(sys.dimension);
  F.alloc(sys.dimension);
  TMAT.alloc(sys.dimension, tsteps);
  TMAT.zeros();
  
  x[0] = 1.06512037300928485983;
  x[1] = 0.51819992367912581788;
  x[2] = 0.85058383242985102779;
  x[3] = -0.66197500054304025952;
  x[4] = -0.01000000000000000;
  x[5] = 0.73618306350367335167;
  x[6] = 0.77067836274882195458;
  x[7] = 0.69832289180288620312;
  x[8] = 0.46185376441989828278;
  x[9] = 1.01531727676021699125;

  x[sys.pnet] = 1.04;
  x[sys.pnet + 1] = 0.0; 
  x[sys.pnet + 2] = 1.01613;
  x[sys.pnet + 3] = -0.05803568828731545;

 
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
  pTensor4<double> T(dim, dim, dim, dim);
  pTensor3<double> H(dim, dim, dim);
  pMatrix<double>  J(dim, dim);

  for (size_t i = 0; i < sys.dimension; ++i) 
    cv0[i][i] = 0.0000001;
  
  cv0[4][4] = 0.000001;

  for (size_t i = 0; i < tsteps; ++i) {
    // Save mean in trajectory matrix
    for (size_t j = 0; j < sys.dimension; ++j) {
      TMAT.set(j, i, x[j]);
    }

    // Propagate
    std::cout << "Step: " << i << ". Time: " << sys.deltat * i
      << "." << std::endl;
    int degree = 3;
    propagateAD(x, cv0, sys, J, H, T, drivers, degree);
  }
  
  std::cout << global_prof << std::endl;
  double compare, diff;
  std::ifstream infile;
  infile.open("solution.txt");
  infile >> compare;
  diff = compare - cv0.norm();
  std::cout << std::setprecision(16) << compare << std::endl;
  std::cout << std::setprecision(16) << cv0.norm() << std::endl;
  std::cout << std::setprecision(16) << compare - cv0.norm() << std::endl;
  infile.close();
  if(diff > 1e-15) {
    cerr << "ERROR: Covariance not correct." << endl;
    exit(1);
  }

  return 0;
}

