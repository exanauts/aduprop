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

  sys.ic(xold);
  ad drivers(sys);


  // jacobian test
  if (test_jac) {
    drivers.jactest(xold);
    exit(0);
  }

  std::cout << "Computing sensitivities at x:" << std::endl;
  std::cout << "---------" << std::endl;
  std::cout << xold << endl;
  x = xold;


  // Tensors
  pMatrix<double>  J(dim, dim);
  pMatrix<double>  J_fd(dim, dim);

  pTensor3<double> H(dim, dim, dim);
  pTensor3<double> H_fd(dim, dim, dim);

  pTensor4<double> T(dim, dim, dim, dim);
  pTensor4<double> T_fd(dim, dim, dim, dim);


  if (tensor1) {
    J.zeros();
    drivers.t1s_driver(x, J);
    std::cout << "Jacobian using 1st order AD" << endl;
    std::cout << "-----------------" << std::endl;
    std::cout << J << std::endl;

    J_fd.zeros();
    drivers.fdJ_driver(x, J_fd);
    std::cout << "Jacobian using FD" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "J" << std::endl << J << std::endl;
    double diff = (J - J_fd).norm();
    cout << "norm(J - J_fd)" << endl << diff << endl;
    if(diff > 1e-4) {
      cerr << "ERROR: HC and AD Jacobian differ." << endl;
      exit(1);
    }
  }

  if (tensor2) {
    J.zeros();
    H.zeros();
    drivers.t2s_t1s_driver(x, J, H);
    std::cout << "Hessian and Jacobian using 2nd order AD" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "J" << std::endl << J << std::endl;
    std::cout << "H" << std::endl << H << std::endl;

    H_fd.zeros();
    std::cout << "Hessian using FD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.fdH_driver(x, H_fd);
    std::cout << "H_fd" << std::endl << H_fd << std::endl;
    
    double diff = (H - H_fd).norm();
    cout << "norm(H - H_fd)" << endl << diff << endl;
    if(diff > 1e-4) {
      cerr << "ERROR: HC and AD Hessian differ." << endl;
      exit(1);
    }
  }

  if (tensor3) {
    J.zeros();
    H.zeros();
    T.zeros();
    std::cout << "Tensor, Hessian and Jacobian using 3rd order AD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.t3s_t2s_t1s_driver(x, J, H, T);
    std::cout << "J" << std::endl << J << std::endl;
    std::cout << "H" << std::endl << H << std::endl;
    std::cout << "T" << std::endl << T << std::endl;
    std::cout << std::endl;

    T_fd.zeros();
    std::cout << "Tensor using FD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.fdT_driver(x, T_fd);
    std::cout << "T_fd" << std::endl << T_fd << std::endl;
    double diff = (T - T_fd).norm();
    cout << "norm(T - T_fd)" << endl << diff << endl;
    if(diff > 1e-5) {
      cerr << "ERROR: HC and AD tensor differ." << endl;
      exit(1);
    }
  }


  return 0;
}

