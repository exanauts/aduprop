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
  System sys;
  size_t dim = sys.dim();

  pVector<double> xold(dim);
  pVector<double> x(dim);

  ad drivers(sys);


  return 0;
}

