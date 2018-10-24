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
#include "parallel.hpp"
#include "mc.hpp"

int main(int argc, char* argv[]) {
  // Options parser


  // activate timer in propagateAD
  
  global_prof.activate("propagateAD");
  global_prof.activate("propagateCOV");
  global_prof.activate("propagateH");
  global_prof.activate("propagateMU");
  global_prof.activate("reduction");
  cxxopts::Options options("UQ Power", "Perform UQ on power system with AD");

  bool do_mc = false;
  bool test_jac = false;
  bool test_tensor1  = false;
  bool test_tensor2  = false;
  bool test_tensor3  = false;
  bool tensor1  = false;
  bool tensor2  = false;
  bool tensor3  = false;

  options.add_options()
    ("test_jac", "Test inner jacobian", cxxopts::value<bool>(test_jac))
    ("mc", "Monte Carlo", cxxopts::value<bool>(do_mc))
    ("test_tensor1", "Compute AD and HC jacobian for 1 ts", cxxopts::value<bool>(test_tensor1))
    ("test_tensor2", "Compute AD and HC hessian for 1 ts", cxxopts::value<bool>(test_tensor2))
    ("test_tensor3", "Compute AD and HC tensor of order 3 for 1 ts",
     cxxopts::value<bool>(test_tensor3))
    ("tensor1", "Compute AD and HC jacobian", cxxopts::value<bool>(tensor1))
    ("tensor2", "Compute AD and HC hessian", cxxopts::value<bool>(tensor2))
    ("tensor3", "Compute AD and HC tensor of order 3",
     cxxopts::value<bool>(tensor3));

  auto result = options.parse(argc, argv);
  // Variable declaration
  // pVector<double> xold, x, F;
  // pVector<double> x0;
  // pMatrix<double> TMAT;
  System sys;
  if (argc < 3)
  {
    std::cout << "plorenz [options] [dimension N] [forcing F] [timesteps]" << std::endl;
    std::cout << "Example: mpiexec -n 2 ./plorenz --tensor1 7 4.4 40000" << std::endl;
    std::cout << "--tensorX: Use x-th (= 1,2 or 3) order derivatives" << std::endl;
    return 1;
  }
  sys.dimension = atoi(argv[1]);
  sys.F = atof(argv[2]);
  size_t tsteps = atoi(argv[3]);
  size_t dim = sys.dim();

  pVector<double> xold(dim);
  pVector<double> x(dim);
  // F.alloc(sys.dimension);

  // Read initial conditions from external file. Check for dimension
  // consistency.
  
  sys.ic(xold);
  ad drivers(sys);
  dim = sys.dimension;
  size_t chunk = paduprop_getend(dim) - paduprop_getstart(dim);
  size_t start = paduprop_getstart(dim);
  size_t end = paduprop_getend(dim);
  x = xold;
  // jacobian test
  if (test_jac)
  {
    drivers.jactest(xold);
    return 0;
  }
  if (test_tensor1) {
    pMatrix<double> J(dim, dim);
    pMatrix<double> J_fd(dim, dim);
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
      return 1;
    }
    return 0;
  }
  if (test_tensor2) {
    pMatrix<double> J(dim, dim);
    pTensor3<double> H(dim, dim, chunk);
    pTensor3<double> H_fd(dim, dim, chunk);
    J.zeros();
    H.zeros();
    drivers.t2s_t1s_driver(x, J, H, start, end);
    std::cout << "Hessian and Jacobian using 2nd order AD" << std::endl;
    std::cout << "-----------------" << std::endl;
    // std::cout << "H" << std::endl << H << std::endl;

    H_fd.zeros();
    drivers.fdH_driver(x, H_fd, start, end);
    
    double diff = (H - H_fd).norm();
    diff = diff*diff;
    paduprop_sum(diff);
    cout << "norm(H - H_fd)" << endl << sqrt(diff) << endl;
    if(diff > 1e-4) {
      cerr << "ERROR: HC and AD Hessian differ." << endl;
      return 1;
    }
    return 0;
  }
  if (test_tensor3) {
    pMatrix<double>  J(dim, dim);
    pTensor3<double> H(dim, dim, dim);
    pTensor4<double> T(dim, dim, dim, chunk);
    pTensor4<double> T_fd(dim, dim, dim, chunk);
    J.zeros();
    H.zeros();
    T.zeros();
    std::cout << "Tensor, Hessian and Jacobian using 3rd order AD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.t3s_t2s_t1s_driver(x, J, H, T, start, end);
    // std::cout << "T" << std::endl << T << std::endl;
    // std::cout << std::endl;

    T_fd.zeros();
    drivers.fdT_driver(x, T_fd, start, end);
    double diff = (T - T_fd).norm();
    diff = diff*diff;
    paduprop_sum(diff);
    cout << "norm(T - T_fd)" << endl << sqrt(diff) << endl;
    if(diff > 1e-5) {
      cerr << "ERROR: HC and AD tensor differ." << endl;
      return 1;
    }
    return 0;
  }
  // Co-variance
  pMatrix<double> cv0(sys.dimension, sys.dimension);
  cv0.zeros();
  // Strings
  std::string fcov;

  // Where do we put this???
  // We should only alocate if we're using all of these.
  pTensor4<double> T;
  pTensor3<double> H;
  pMatrix<double>  J(dim, dim);
  if (tensor2)
  {
    H=pTensor3<double>(dim, dim, chunk);
  }
  if (tensor3)
  {
    T=pTensor4<double>(dim, dim, dim, chunk);
    H=pTensor3<double>(dim, dim, dim);
  }

  for (size_t i = 0; i < sys.dimension; ++i) 
    cv0[i][i] = 0.0000001;
  
  //cv0[sys.dim()][sys.dimension] = 0.000001;
  if (do_mc)
  {
    size_t samples=1000;
    mc mc_integration(samples, sys); 
    mc_integration.integrate(x, cv0, tsteps);
    return 0;
  }

  std::ofstream xfile;
  std::ofstream covfile;
  if (paduprop_getrank() == 0)
  {
    xfile.open("datas");
    if (tensor2)
    {
      covfile.open("data_cov2");
    }
    else
    {
      if (tensor3)
      {
        covfile.open("data_cov3");
      }
      else
      {
        covfile.open("data_cov1");
      }
    }
  }
  for (size_t i = 0; i < tsteps; ++i) {
    // Save mean in trajectory matrix
    // for (size_t j = 0; j < sys.dimension; ++j) {
    //   TMAT.set(j, i, x[j]);
    // }

    // Propagate
    if (paduprop_getrank() == 0)
    {
      std::cout << "Step: " << i << "." << std::endl;
    }
    int degree = 1;
    if(tensor1) degree = 1;
    if(tensor2) degree = 2;
    if(tensor3) degree = 3;
    propagateAD(x, cv0, sys, J, H, T, drivers, degree);
    if (paduprop_getrank() == 0)
    {
      // xfile << "X: " << i * sys.h << " " << x << std::endl;
      // covfile << "COV: " << i * sys.h << " ";
      for (size_t j = 0; j < cv0.nrows(); j++)
      {
        xfile << x[j] << " ";
        covfile << cv0[j][j] << " ";
      }
      xfile << std::endl;
      covfile << std::endl;
    }
      // xfile << "X: " << i * sys.h << " " << x << std::endl;
      // covfile << "COV: " << i * sys.h << " ";
      for (size_t j = 0; j < cv0.nrows(); j++)
      {
        std::cout << x[j] << " ";
      }
      std::cout << std::endl;
      for (size_t j = 0; j < cv0.nrows(); j++)
      {
        std::cout << cv0[j][j] << " ";
      }
      std::cout << std::endl;
  }
  // cout << cv0;
  if(paduprop_getrank() == 0) {  
    xfile.close();
    covfile.close();
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
}

