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

void test(int argc, char* argv[], pVector<double> xold, System& sys,
    ad& drivers) {
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
  size_t dim = sys.dim();
  
  pVector<double> x(dim);

  // jacobian test
  if (test_jac) {
    drivers.jactest(xold);
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
  }

  if (tensor2) {
    J.zeros();
    H.zeros();
    drivers.t2s_t1s_driver(x, J, H);
    std::cout << "Hessian and Jacobian using 2nd order AD" << std::endl;
    std::cout << "-----------------" << std::endl;
    std::cout << "J" << std::endl << J << std::endl;
    std::cout << "H" << std::endl << H << std::endl;

    J.zeros();
    H.zeros();
    std::cout << "Hessian using FD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.fdH_driver(x, H);
    std::cout << H << std::endl;
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

    J.zeros();
    H.zeros();
    T.zeros();
    std::cout << "Tensor using FD" << std::endl;
    std::cout << "-----------------" << std::endl;
    drivers.fdT_driver(x, T);
    std::cout << T << std::endl;
  }
}

// TODO(Michel, Adrian): names for system objects and ad drivers seem to me a tad
// confusing. Can we think how to name this?

/*!
   \brief "Propagates the mean and the covariance matrices one step."
   \param m0 "Mean"
   \param cv0 "Covariance matrix"
   \param sys "System data"
   \param drivers "AD drivers"
*/
void propagateAD(pVector<double>& m0, pMatrix<double>& cv0, System& sys,
    ad& drivers) {
  
  size_t dim = sys.dim();
  pMatrix<double>  J(dim, dim);
  pMatrix<double>  cv_temp(dim, dim);
  

  // Before
  //std::cout << m0 << std::endl;
  //std::cout << cv0 << std::endl;
  
  // Obtain tensors
  J.zeros();
  cv_temp.zeros();
  drivers.t1s_driver(m0, J);
  
  // Propagate mean
  drivers.integrate(m0);


  // Propagate covariance

  for (size_t pn = 0; pn < dim; ++pn) {
    for (size_t pm = 0; pm < dim; ++pm) {
      for (size_t i = 0; i < dim; ++i) {
        for (size_t j = 0; j < dim; ++j) {
          cv_temp[pn][pm] += 0.5*((J[pn][j]*J[pm][i] +
                J[pm][j]*J[pn][i])*cv0[i][j]);
        }
      }
    }
  }

  cv0 = cv_temp;

  // After
  //std::cout << m0 << std::endl;
  //std::cout << cv0 << std::endl;

}



int main(int argc, char* argv[]) {
  
  // problem definition
  System sys;
  size_t dim = sys.dim();

  // Mean
  pVector<double> m0(dim);
  m0[0] = 1.06512037300928485983;
  m0[1] = 0.51819992367912581788;
  m0[2] = 0.85058383242985102779;
  m0[3] = -0.66197500054304025952;
  m0[4] = -0.05000000000000000278;
  m0[5] = 0.73618306350367335167;
  m0[6] = 0.77067836274882195458;
  m0[7] = 0.69832289180288620312;
  m0[8] = 0.46185376441989828278;
  m0[9] = 1.01531727676021699125;
  m0[10] = 1.0400000000000000000;
  m0[11] = 0.00000000000000000000;
  
  // Co-variance
  pMatrix<double> cv0(dim, dim);
  cv0.zeros();
  for (size_t i = 0; i < dim; ++i) {
    cv0[i][i] = 0.0001;
  }


  ad drivers(sys);
  
  //test(argc, argv, m0, sys, drivers);

  for (size_t i = 0; i < 30; ++i) {
    std::cout << "Step: " << i << std::endl;
    propagateAD(m0, cv0, sys, drivers);
  }


  return 0;
}

