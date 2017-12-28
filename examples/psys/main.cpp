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
// AD library used
// Rudimentary linear solver interface
#include "ad.hpp"
#include "cxxopts.hpp"

int main(int argc, char* argv[]) {
  // Options parser

  cxxopts::Options options("UQ Power", "Perform UQ on power system with AD");

  bool test_jac = false;
  bool tensor1  = false;
  bool tensor2  = false;
  bool tensor3  = false;

  options.add_options()
    ("test_jac", "Test inner jacobian", cxxopts::value<bool>(test_jac))
    ("tensor1", "Compute AD and HC jacobian", cxxopts::value<bool>(tensor1))
    ("tensor2", "Compute AD and HC hessian", cxxopts::value<bool>(tensor2))
    ("tensor3", "Compute AD and HC tensor of order 3",
     cxxopts::value<bool>(tensor3));

  auto result = options.parse(argc, argv);

  // problem definition
  int nbuses, nbranches, ngen, nload;

  nbuses = 2;
  nbranches = 1;
  ngen = 1;
  nload = 1;

  System sys(nbuses, nbranches, ngen, nload);

  sys.branches[0].set(0, 1, 0.0001, 0.0576);
  sys.loads[0].set(1, 1.0648453, 0.38835684);
  sys.gens[0].set(0, 1.575, 1.512, 0.29, 0.39, 0.1733,
      0.1733, 0.0787, 3.38, 6.09, 1.0, 0.05, 0.15);
  sys.gens[0].e_fd = 2.36980307364616349375;
  sys.gens[0].p_m = 1.06496000000000012875;

  sys.build_ybus();
  //std::cout << *sys.ybus << std::endl;

  pVector<double> xold(sys.dimension);
  pVector<double> x(sys.dimension);
  pVector<double> F(sys.dimension);

  x[0] = 1.04;
  x[1] = 0.0; 
  x[2] = 1.01613;
  x[3] = -0.05803568828731545;

  xold = x;

  sys.residual_beuler<double>(x, xold, F);
  std::cout << F << std::endl;

  //ad drivers(sys);


  return 0;
}

