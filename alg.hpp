#ifndef ADUPROP_ALG_HPP_
#define ADUPROP_ALG_HPP_
/*!
   \file "alg.hpp"
   \brief "Linear algebra module"
   \author "Adrian Maldonado and Michel Schanen"
   \date 21/11/2017
*/

#include <codi.hpp>
#include <cassert>
#include "tensor.hpp"

namespace alg {

// First, second and third order AD datatype
typedef codi::RealForwardGen<double> t1s;
typedef codi::RealForwardGen<t1s> t2s;
typedef codi::RealForwardGen<t2s> t3s;

// Passive vector definition

template <typename> class pVector;
template <typename T> std::ostream& operator<< (std::ostream& os,
    pVector<T> &m);

template <typename T> class pVector {
 public:
  pVector();
  explicit pVector(const size_t nvals);
  ~pVector();

  void set(const size_t i, const T val);
  T get(const size_t i) const;
  size_t dim() const;
  T* get_datap() const;
  void zeros();

#ifdef HDF5
  // IO routines
  void to_hdf5(const std::string filename);
  void from_hdf5(const std::string filename);
#endif

  void display();
  T& operator[] (int i) {
    return data[i];
  }

  const T& operator[] (int i) const {
    return data[i];
  }

  friend std::ostream& operator<< <> (std::ostream&, pVector<T>&);

  pVector(const pVector &other) {
    if (other.n != n) {
      if (data != NULL) delete [] data;
      n = other.n;
      data = new T[n];
    }
    for (size_t i = 0; i < n; ++i) data[i] = other.data[i];
  }

  pVector& operator=(const pVector &other) {
    if (other.n != n) {
      if (data != NULL) delete [] data;
      n = other.n;
      data = new T[n];
    }
    for (size_t i = 0; i < n ; ++i) data[i] = other.data[i];
    return *this;
  }

  pVector operator+(const pVector& b) {
    assert(this->n == b.n);
    pVector<T> vec(b.n);
    for (size_t i = 0; i < b.n; ++i) {
      vec.data[i] = this->data[i] + b.data[i];
    }
    return vec;
  }

  pVector operator-(const pVector& b) {
    assert(this->n == b.n);
    pVector<T> vec(b.n);
    for (size_t i = 0; i < b.n; ++i) {
      vec.data[i] = this->data[i] - b.data[i];
    }
    return vec;
  }

  T operator*(const pVector& b) {
    assert(this->n == b.n);
    T res = 0.0;
    for (size_t i = 0; i < b.n; ++i) {
      res += this->data[i] * b.data[i];
    }
    return res;
  }
  T norm() {
    return sqrt((*this) * (*this));
  }

 private:
  size_t n;
  T* data = NULL;
};

template <typename T> inline size_t pVector<T>::dim() const {
  return n;
}

template <typename T> inline T* pVector<T>::get_datap() const {
  return data;
}

template <typename T> inline void pVector<T>::zeros() {
  for (size_t i = 0; i < n; ++i) {
    data[i] = 0.0;
  }
}

template <typename T> std::ostream& operator<< (std::ostream& os,
    pVector<T> &v) {
  os << "[ ";
  for (size_t i = 0; i < v.n; ++i) {
    if (i) os << ", ";
    os << v.data[i];
  }
  os << " ]";
  std::cout << std::endl;
  return os;
}

// Passive matrix definition
template <typename> class pMatrix;
template <typename T> std::ostream& operator<< (std::ostream& os,
    pMatrix<T> &m);

template <typename T> class pMatrix {
 public:
  pMatrix();
  explicit pMatrix(const size_t nrows, const size_t ncols);
  ~pMatrix();

  void set(const size_t i, const size_t j, const T val);
  void set_col(const size_t j, const T *vals);
  T& get(const size_t i, const size_t j);
  T* get_datap() const;
  void zeros();
  size_t nrows() const;
  size_t ncols() const;
  void display();

  class row {
   public:
    T *ptr;
    size_t &rows;
    row(T *ptr_, size_t &rows_) : ptr(ptr_), rows(rows_) {};
    T& operator[] (int i) {
      return *(ptr+i*rows);
    }
  };

  row operator[] (int i) {
    return row(data+i, rows);
  }

  friend std::ostream& operator<< <> (std::ostream&, pMatrix<T>&);

  pMatrix(const pMatrix &other) {
    if (other.cols != cols || other.rows != rows) {
      if (data != NULL) delete [] data;
      rows = other.rows;
      cols = other.cols;
      data = new T[rows*cols];
    }
    for (size_t i = 0; i < rows*cols; ++i) data[i] = other.data[i];
  }
  pMatrix& operator=(const pMatrix &other) {
    if (other.cols != cols || other.rows != rows) {
      if (data != NULL) delete [] data;
      rows = other.rows;
      cols = other.cols;
      data = new T[rows*cols];
    }
    for (size_t i = 0; i < rows*cols; ++i) data[i] = other.data[i];
    return *this;
  }

  T norm() {
    T res = 0;
    for (size_t i = 0 ; i < rows; ++i) {
      for (size_t j = 0 ; j < cols; ++j) {
        res+=(data[i*cols + j])*(data[i*cols + j]);
      }
    }
    return sqrt(res);
  }
  pVector<T> operator*(const pVector<T>& b) {
    assert(this->cols == b.dim());
    pVector<T> res(rows);
    for (size_t i = 0; i < rows; ++i) res[i] = 0;
    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        res[i] += data[j*cols + i] * b[j];
      }
    }
    return res;
  }
  
 private:
  size_t rows, cols;
  T* data;
};

template <typename T> inline void pMatrix<T>::zeros() {
  for (size_t i = 0; i < cols*rows; ++i) {
    data[i] = 0.0;
  }
}

template <typename T> std::ostream& operator<< (std::ostream& os,
    pMatrix<T> &m) {
  for (size_t i = 0; i < m.rows; ++i) {
    os << "[ ";
    for (size_t j = 0; j < m.cols; ++j) {
      if (j) os << ", ";
      os << m.data[j*m.cols + i];
    }
    os << " ],\n";
  }
  return os;
}

template <typename T> inline size_t pMatrix<T>::nrows() const {
  return rows;
}

template <typename T> inline size_t pMatrix<T>::ncols() const {
  return cols;
}

template <typename T> inline T* pMatrix<T>::get_datap() const {
  return data;
}

template <typename T> inline T& pMatrix<T>::get(const size_t i,
    const size_t j) {
  return data[j*rows + i];
}

// Function declarations

void decmatmul(const pMatrix<double> &A, const pVector<double> &x,
  pVector<double> &y);

void LUsolve(pMatrix<double> &A, pVector<double> &b);

}  // namespace alg

#endif  // ADUPROP_ALG_HPP_
