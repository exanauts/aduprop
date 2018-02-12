#ifndef ADUPROP_TENSOR_HPP_
#define ADUPROP_TENSOR_HPP_
/*!
   \file "alg.hpp"
   \brief "Linear algebra module"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <codi.hpp>
#include <cassert>

namespace alg {

// First, second and third order AD datatype
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

// Tensor for 3 and 4 dimensions
template <typename T> class pTensor3;
template <typename T> std::ostream& operator<< (std::ostream& os, pTensor3<T> &m);

template <typename T> class pTensor3 {
 public:

  class cd1 {
  public:
    class cd2 {
    public:
      T *ptr;
      cd2(T *ptr_) : ptr(ptr_) {}; 
      T& operator[] (int i) {
        return *(ptr+i);
      }
      
    };
    T *ptr;
    size_t &d1;
    cd1(T *ptr_, size_t &d1_) : ptr(ptr_), d1(d1_) {}; 
    cd2 operator[] (int i) {
      return cd2(ptr+i*d1);
    }
  
  };
  cd1 operator[] (int i) {
    return cd1(data+i*d2*d3, d3);
  }
  friend std::ostream& operator<< <> ( std::ostream&, pTensor3<T>& );  
  
  pTensor3( const pTensor3 &other ) {
    if(d1 != other.d1 || d2 != other.d2 || d3 != other.d3) {
      if(data != NULL) delete [] data;
      d1=other.d1;
      d2=other.d2;
      d3=other.d3;
      data = new T[d1*d2*d3];
    }
    for(size_t i = 0; i < d1*d2*d3 ; ++i) data[i] = other.data[i];
  }
  pTensor3<T> operator+(const pTensor3<T>& b) {
    assert(this->d1 == b.d1);
    assert(this->d2 == b.d2);
    assert(this->d3 == b.d3);
    pTensor3<T> res(d1,d2,d3);
    T *pb = b.get_datap();
    T *pres = res.get_datap();
    for (size_t i = 0; i < d1*d2*d3; ++i) {
        pres[i] = data[i] + pb[i];
    }
    return res;
  }
  
  pTensor3<T> operator-(const pTensor3<T>& b) {
    assert(this->d1 == b.d1);
    assert(this->d2 == b.d2);
    assert(this->d3 == b.d3);
    pTensor3<T> res(d1,d2,d3);
    T *pb = b.get_datap();
    T *pres = res.get_datap();
    for (size_t i = 0; i < d1*d2*d3; ++i) {
        pres[i] = data[i] - pb[i];
    }
    return res;
  }
  
  pTensor3& operator=( const pTensor3 &other ) {
    if(d1 != other.d1 || d2 != other.d2 || d3 != other.d3) {
      if(data != NULL) delete [] data;
      d1=other.d1;
      d2=other.d2;
      d3=other.d3;
      data = new T[d1*d2*d3];
    }
    for(size_t i = 0; i < d1*d2*d3 ; ++i) data[i] = other.data[i];
    return *this;
  }
  
  T norm() {
    T res = 0;
    for(size_t i = 0 ; i < d1 ; ++i) {
      for(size_t j = 0 ; j < d2 ; ++j) {
        for(size_t k = 0 ; k < d3 ; ++k) {
          res+=(data[i*d2*d3 + j*d3 + k] * data[i*d2*d3 + j*d3 + k]);
        }
      }
    }
    return sqrt(res);
  }
  
  
  void zeros() {
    for (size_t i = 0; i < d3*d2*d1; ++i) {
      data[i] = 0.0;
    }
  }
  
  size_t get_d1() const {
    return d1;
  }
  
  size_t get_d2() const {
    return d2;
  }
  
  size_t get_d3() const {
    return d3;
  }
  
  T* get_datap() const {
    return data;
  }
  
  pTensor3() {
    data = NULL;
    d3 = 0;
    d2 = 0;
    d1 = 0;
  }
  
  pTensor3(const size_t d1_, const size_t d2_, const size_t d3_) {
    data = new T[d3_*d2_*d1_];
    d3 = d3_;
    d2 = d2_;
    d1 = d1_;
  }
  
  ~pTensor3() {
    if (data != NULL) {
      delete [] data;
      data = NULL;
    }
  }
  
private:
  size_t d1, d2, d3;
  T* data;
};


template <typename T> std::ostream& operator<< (std::ostream& os, alg::pTensor3<T> &m) {
  for (size_t i = 0; i < m.d1; ++i) {
    for (size_t j = 0; j < m.d2; ++j) {
      for (size_t k = 0; k < m.d3; ++k) {
        os << m[i][j][k] << " ";
      }
    }
    os << "\n";
  }
  return os;
}

template <typename T> class pTensor4;
template <typename T> std::ostream& operator<< (std::ostream& os, pTensor4<T> &m);

template <typename T> class pTensor4 {
 public:

  class cd1 {
  public:
    class cd2 {
    public:
      class cd3 {
      public:
        T *ptr;
        cd3(T *ptr_) : ptr(ptr_) {}; 
        T& operator[] (int i) {
          return *(ptr+i);
        }
        
      };
      T *ptr;
      size_t &d1;
      cd2(T *ptr_, size_t &d1_) : ptr(ptr_), d1(d1_) {}; 
      cd3 operator[] (int i) {
        return cd3(ptr+i*d1);
      }
      
    };
    T *ptr;
    size_t &d1;
    size_t &d2;
    cd1(T *ptr_, size_t &d1_, size_t &d2_) : ptr(ptr_), d1(d1_), d2(d2_) {}; 
    cd2 operator[] (int i) {
      return cd2(ptr+i*d1*d2, d2);
    }
  
  };
  cd1 operator[] (int i) {
    return cd1(data+i*d2*d3*d4, d3, d4);
  }
  friend std::ostream& operator<< <> ( std::ostream&, pTensor4<T>& );  
  
  pTensor4( const pTensor4 &other ) {
    if(d1 != other.d1 || d2 != other.d2 || d3 != other.d3 || d4 != other.d4) {
      if(data != NULL) delete [] data;
      d1=other.d1;
      d2=other.d2;
      d3=other.d3;
      d3=other.d4;
      data = new T[d1*d2*d3*d4];
    }
    for(size_t i = 0; i < d1*d2*d3*d4 ; ++i) data[i] = other.data[i];
  }
  
  pTensor4<T> operator+(const pTensor4<T>& b) {
    assert(this->d1 == b.d1);
    assert(this->d2 == b.d2);
    assert(this->d3 == b.d3);
    assert(this->d4 == b.d4);
    pTensor4<T> res(d1,d2,d3,d4);
    T *pb = b.get_datap();
    T *pres = res.get_datap();
    for (size_t i = 0; i < d1*d2*d3*d4; ++i) {
        pres[i] = data[i] + pb[i];
    }
    return res;
  }
  
  pTensor4<T> operator-(const pTensor4<T>& b) {
    assert(this->d1 == b.d1);
    assert(this->d2 == b.d2);
    assert(this->d3 == b.d3);
    assert(this->d4 == b.d4);
    pTensor4<T> res(d1,d2,d3,d4);
    T *pb = b.get_datap();
    T *pres = res.get_datap();
    for (size_t i = 0; i < d1*d2*d3*d4; ++i) {
        pres[i] = data[i] - pb[i];
    }
    return res;
  }
  
  pTensor4& operator=( const pTensor4 &other ) {
    if(d1 != other.d1 || d2 != other.d2 || d3 != other.d3 || d4 != other.d4) {
      if(data != NULL) delete [] data;
      d1=other.d1;
      d2=other.d2;
      d3=other.d3;
      d4=other.d4;
      data = new T[d1*d2*d3*d4];
    }
    for(size_t i = 0; i < d1*d2*d3*d4 ; ++i) data[i] = other.data[i];
    return *this;
  }
  
  T norm() {
    T res = 0;
    for(size_t i = 0 ; i < d1 ; ++i) {
      for(size_t j = 0 ; j < d2 ; ++j) {
        for(size_t k = 0 ; k < d3 ; ++k) {
          for(size_t l = 0 ; l < d4 ; ++l) {
            res+=(data[i*d2*d3*d4 + j*d3*d4 + k*d4 + l] * data[i*d2*d3*d4 + j*d3*d4 + k*d4 + l]);
          }
        }
      }
    }
    return sqrt(res);
  }
  
  
  void zeros() {
    for (size_t i = 0; i < d4*d3*d2*d1; ++i) {
      data[i] = 0.0;
    }
  }
  
  size_t get_d1() const {
    return d1;
  }
  
  size_t get_d2() const {
    return d2;
  }
  
  size_t get_d3() const {
    return d3;
  }
  
  size_t get_d4() const {
    return d4;
  }
  
  T* get_datap() const {
    return data;
  }
  
  pTensor4() {
    data = NULL;
    d4 = 0;
    d3 = 0;
    d2 = 0;
    d1 = 0;
  }
  
  pTensor4(const size_t d1_, const size_t d2_, const size_t d3_, const size_t d4_) {
    data = new T[d4_*d3_*d2_*d1_];
    d4 = d4_;
    d3 = d3_;
    d2 = d2_;
    d1 = d1_;
  }
  
  ~pTensor4() {
    if (data != NULL) {
      delete [] data;
      data = NULL;
    }
  }
  
private:
  size_t d1, d2, d3, d4;
  T* data;
};


template <typename T> std::ostream& operator<< (std::ostream& os, alg::pTensor4<T> &m) {
  for (size_t i = 0; i < m.d1; ++i) {
    for (size_t j = 0; j < m.d2; ++j) {
      for (size_t k = 0; k < m.d3; ++k) {
        for (size_t l = 0; l < m.d4; ++l) {
          os << m[i][j][k][l] << " ";
        }
      }
    }
    os << "\n";
  }
  return os;
}

} // end of namespace
#endif  // ADUPROP_TENSOR_HPP_
