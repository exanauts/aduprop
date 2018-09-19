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

  bool external_init = false;
  bool propagate_moments = false;

  bool TWO_BUS = false;
  int tsteps = 200;
  
  options.add_options()
    ("x0", "Read initial conditions from x0.hdf5", cxxopts::value<bool>(external_init))
    ("m,moments", "Propagate moments with AD", cxxopts::value<bool>(propagate_moments))
    ("o,output", "Output trajectory (integration)", cxxopts::value<std::string>())
    ("t,timesteps", "Number of time steps", cxxopts::value<int>(tsteps));

  auto result = options.parse(argc, argv);

  // activate timer in propagateAD
  
  global_prof.activate("propagateAD");

  // Variable declaration
  int nbuses, nbranches, ngen, nload;
  //int tsteps = 200;
  size_t dim;
  pVector<double> xold, x, F;
  pVector<double> x0;
  pMatrix<double> TMAT;
  System sys;


  if (TWO_BUS) {
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

  } else {

    nbuses = 9;
    nbranches = 9;
    ngen = 3;
    nload = 3;

    sys.init(nbuses, nbranches, ngen, nload);
    
    sys.branches[0].set(0, 3, 0.0000, 0.0576, 0.000);
    sys.branches[1].set(1, 6, 0.0000, 0.0625, 0.000);
    sys.branches[2].set(2, 8, 0.0000, 0.0586, 0.000);
    sys.branches[3].set(3, 4, 0.0100, 0.0850, 0.176);
    sys.branches[4].set(3, 5, 0.0170, 0.0920, 0.158);
    sys.branches[5].set(4, 6, 0.0320, 0.1610, 0.306);
    sys.branches[6].set(5, 8, 0.0390, 0.1700, 0.358);
    sys.branches[7].set(6, 7, 0.0085, 0.0720, 0.149);
    sys.branches[8].set(7, 8, 0.0119, 0.1008, 0.209);

    // sys.loads[0].set(4, 1.383061, 0.553224);
    sys.loads[0].set(4, 1.2500001068394062, 0.49999968121877764);

    //sys.loads[1].set(5, 0.964050, 0.321350);
    sys.loads[1].set(5, 0.900000218680605,  0.300000072893535);
    sys.loads[2].set(7, 0.9999996236911225, 0.3500004914694625);

    // Generator 0
    sys.gens[0].set(0, 1.575, 1.512, 0.29, 0.39, 0.1733,
        0.1733, 0.0787, 3.38, 6.09, 1.0, 0.05, 0.15);
    sys.gens[0].e_fd = 2.52242991977;
    sys.gens[0].p_m = 0.7234;

    sys.gens[1].set(1, 1.575, 1.512, 0.29, 0.39, 0.1733,
        0.1733, 0.0787, 3.38, 6.09, 1.0, 0.05, 0.15);
    sys.gens[1].e_fd = 3.11103600578;
    sys.gens[1].p_m = 1.63;

    sys.gens[2].set(2, 1.575, 1.512, 0.29, 0.39, 0.1733,
        0.1733, 0.0787, 3.38, 6.09, 1.0, 0.05, 0.15);
    sys.gens[2].e_fd = 2.07873175984;
    sys.gens[2].p_m = 0.85;

    sys.build_ybus();
    x.alloc(sys.dimension);
    F.alloc(sys.dimension);

    TMAT.alloc(sys.dimension, tsteps);
    TMAT.zeros();

    // Generator 0 state variables

    x[0] = 1.23144118e+00;
    x[1] = 3.30056857e-01;
    x[2] = 1.01915642e+00;
    x[3] = -4.21631455e-01;
    x[4] = -3.10192730e-25;
    x[4] = 0.0;
    x[5] =  4.41919647e-01;
    x[6] = 0.94008964;
    x[7] =  0.4447825;
    x[8] =  0.29416832;
    x[9] =  1.0046605;
    // Generator 1 state variables
    x[10] = 1.07627648e+00; 
    x[11] = 6.07336946e-01;
    x[12] = 7.41689176e-01;
    x[13] = -7.75843177e-01;
    x[14] = 7.23783036e-25;
    x[15] = 1.09393615e+00;
    x[16] = 0.61707005;
    x[17] = 0.81844337;
    x[18] = 0.54129853;
    x[19] = 1.58347045;

    // Generator 2 state variables
    x[20] = 1.03999032e+00;
    x[21] = 4.70311463e-01;
    x[22] = 8.69184042e-01;
    x[23] = -6.00799839e-01;
    x[24] = 2.03563979e-25;
    x[25] = 7.51809602e-01;

    x[26] = 0.80556618;
    x[27] = 0.63378871;
    x[28] = 0.41917243;
    x[29] = 0.8083591;

    // Bus voltage vector
    x[sys.pnet] = 1.04;
    x[sys.pnet + 1] = 0.0; 
    x[sys.pnet + 2] = 1.02500;
    x[sys.pnet + 3] = (M_PI/180.0)*9.6926;
    x[sys.pnet + 4] = 1.02500;
    x[sys.pnet + 5] = (M_PI/180.0)*4.8812;
    x[sys.pnet + 6] = 0.99574;
    x[sys.pnet + 7] = (M_PI/180.0)*-2.3060;
    x[sys.pnet + 8] = 0.95068;
    x[sys.pnet + 9] = (M_PI/180.0)*-4.1382;
    x[sys.pnet + 10] = 0.96621;
    x[sys.pnet + 11] = (M_PI/180.0)*-3.7372;
    x[sys.pnet + 12] = 0.99740;
    x[sys.pnet + 13] = (M_PI/180.0)*3.9736;
    x[sys.pnet + 14] = 0.97915;
    x[sys.pnet + 15] = (M_PI/180.0)*0.8364;
    x[sys.pnet + 16] = 1.00414;
    x[sys.pnet + 17] = (M_PI/180.0)*2.1073;
  }
 
  // Read initial conditions from external file. Check for dimension
  // consistency.
  if (external_init) {
    x0.from_hdf5("x0.hdf5");
    assert(x0.dim() == x.dim());
    x = x0;
  }

  ad drivers(sys);

  dim = sys.dimension;

  if (!propagate_moments) {
    for (size_t i = 0; i < tsteps; ++i) {
      for (size_t j = 0; j < sys.dimension; ++j) {
        TMAT.set(j, i, x[j]);
      }
      std::cout << "Step: " << i << ". Time: " << sys.deltat * i
        << "." << std::endl;
      drivers.integrate(x);
    }
  } else {
    // Co-variance
    pMatrix<double> cv0(sys.dimension, sys.dimension);
    cv0.zeros();
    // Strings
    char step_str[20];
    std::string fcov;

    // Where do we put this???
    // We should only alocate if we're using all of these.
    size_t chunk = paduprop_getend(dim) - paduprop_getstart(dim);
    std::cout << "chunk: " << chunk << std::endl;
    pTensor4<double> T(dim, dim, dim, chunk);
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
      // Save cov. matrix .... Strings = vodoo
      sprintf(step_str, "%d", (int)i);
      fcov = std::string("output/covSTEP") + step_str + ".hdf5";
      cv0.to_hdf5(fcov.c_str());

      // Propagate
      std::cout << "Step: " << i << ". Time: " << sys.deltat * i
        << "." << std::endl;
      propagateAD(x, cv0, sys, J, H, T, drivers, 3);
    }
  }

  // Output trajectory
  if (result.count("output")) {
    TMAT.to_hdf5(result["output"].as<std::string>());
  } else {
    TMAT.to_hdf5("solution.hdf5");
  }
  
  std::cout << global_prof << std::endl;

  return 0;
}

