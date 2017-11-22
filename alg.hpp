#ifndef ALG_HPP_
#define ALG_HPP_
/*!
   \file "alg.hpp"
   \brief "Linear algebra module"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <codi.hpp>

// First, second and third order AD datatype
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;


// Passive vector definition

class pVector {
  public:
    explicit pVector(size_t nvals);
    ~pVector();
    
    void set(size_t i, double val);
    double get(size_t i);


  private:
    size_t n;
    double* data;
};

// Passive matrix definition

class pMatrix {
  public:
    explicit pMatrix(size_t nrows, size_t ncols);
    ~pMatrix();

    void set(size_t i, size_t j, double val);
    void set_col(size_t j, double *vals);
    double get(size_t i, size_t j);
    void display();

  private:
    size_t rows, cols;
    double* data;
};


// First order active matrix

class t1Matrix {
  public:
};

class t2Matrix {
  public:
};


#endif
