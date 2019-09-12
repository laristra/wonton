/*
 This file is part of the Ristra wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_MATRIX_H_
#define WONTON_MATRIX_H_

#include <cassert>
#include <vector>
#include <array>
#include <iostream>
#include <type_traits>

#include "wonton/support/Vector.h"   // Wonton::Vector
#include "wonton/support/wonton.h"   // Wonton::pow2

/*!
  @file Matrix.h
  @brief Matrix class for Wonton
*/

namespace Wonton {

static std::string ignore_msg="ignore";
static std::string &ignore_msg_ref(ignore_msg);

class Matrix {
 public:
  Matrix() : Rows_(0), Columns_(0) {}
  
  Matrix(int const Rows, int const Columns) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);    // uninitialized
  }
  
  // Initialize to some constant value
  explicit Matrix(int const Rows, int const Columns,
                  double initval) :
      Rows_(Rows), Columns_(Columns) {
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = initval;
  }
  
  // Initialize from vector of vectors - assume that each row vector
  // is of the same size - if its not, then some values will be
  // uninitialized
  
  explicit Matrix(std::vector<std::vector<double>> const& invals) {
    Rows_ = invals.size();
    Columns_ = invals[0].size();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = invals[i][j];
  }

  Matrix(Matrix const& M) : Rows_(M.rows()), Columns_(M.columns()) {
    A_ = M.data();
  }

  Matrix& operator=(Matrix const& M) {
    Rows_ = M.rows();
    Columns_ = M.columns();
    A_.resize(Rows_*Columns_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] = M[i][j];

    is_singular_ = 0;

    return *this;
  }

  /// Add the Matrix @c rhs to this Matrix.
  Matrix& operator+=(const Matrix& rhs) {
    assert((Rows_ == rhs.rows()) && (Columns_ == rhs.columns()));

    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] += rhs[i][j];

    is_singular_ = 0;
  
    return *this;
  }

  /// Subtract the Matrix @c rhs from this Matrix.
  Matrix& operator-=(const Matrix& rhs) {
    assert((Rows_ == rhs.rows()) && (Columns_ == rhs.columns()));

    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        A_[i*Columns_+j] -= rhs[i][j];

    is_singular_ = 0;
  
    return *this;
  }  

  /// Multiply this Matrix by a scalar
  Matrix& operator*=(const double rhs) {
    for (int i = 0; i < A_.size(); ++i)
      A_[i] *= rhs;

    is_singular_ = 0;
  
    return *this;
  }

  // Destructor
  ~Matrix() {}

  //! Number of rows
  int rows() const {return Rows_;}

  //! Number of columns
  int columns() const {return Columns_;}

  //! Return a row of values
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of the matrix
  // Since we don't know what the host code will do with the pointer, 
  // (i.e. change values) we invalidate the singularity indicator.
  
  double * operator[](int const RowIndex) {
    is_singular_ = 0;
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Return a row of values that cannot be modified
  //
  // The main utility of this operator is to facilitate the use
  // of the [][] notation to refer to elements of a const matrix
  
  double const * operator[](int const RowIndex) const {
    return &(A_[RowIndex*Columns_]);
  }
  
  //! Get a transpose
  Matrix transpose() const {
    Matrix AT(Columns_, Rows_);
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < Columns_; ++j)
        AT[j][i] = A_[i*Columns_+j];
    return AT;
  }

  //! Get Inverse - only if its a square matrix

  Matrix inverse() {
    if (Rows_ != Columns_) {
      std::cerr << "Matrix is not rectangular" << std::endl;
      throw std::exception();
    }
    
    Matrix Ainv(Rows_, Rows_, 0.0);

    // Create a temporary matrix with twice as many columns
    Matrix D(Rows_, 2*Rows_);
    
    // Initialize the reduction matrix
    int Rows2 = 2*Rows_;
    for (int i = 0; i < Rows_; i++) {
      for (int j = 0; j < Rows_; j++) {
        D[i][j] = A_[i*Columns_+j];
        D[i][Rows_+j] = 0.0;
      }
      D[i][Rows_+i] = 1.0;
    }
    
    // Do the reduction
    for (int i = 0; i < Rows_; i++) {
      double alpha = D[i][i];
      if (alpha == 0.0) {
        is_singular_ = 2;
        return Ainv;
      }
      
      for (int j = 0; j < Rows2; j++)
        D[i][j] = D[i][j]/alpha;
      
      for (int k = 0; k < Rows_; k++) {
        if ((k - i) == 0)
          continue;
        
        double beta = D[k][i];
        for (int j = 0; j < Rows2; j++)
          D[k][j] = D[k][j] - beta*D[i][j];
      }
    }
    is_singular_ = 1;
    
    // Copy result into output matrix
    for (int i = 0; i < Rows_; i++)
      for (int j = 0; j < Rows_; j++)
        Ainv[i][j] = D[i][j + Rows_];

    return Ainv;
  }
  
  /*!
    @brief  Matrix-Vector multiply with std::vector
    @param[in] X The vector to post-multiply with
    
  */
  
  std::vector<double> operator*(std::vector<double> const& X) const {
    assert(Columns_ == X.size());

    std::vector<double> AX(Rows_);
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Vector multiply with Wonton::Vector
    @param[in] X The vector to post-multiply with
    
  */
  
  template<int D>
  Vector<D> operator*(Vector<D> const& X) const {
    assert(Rows_ == D && Columns_ == D);

    Vector<D> AX;
    for (int i = 0; i < Rows_; ++i) {
      AX[i] = 0.0;
      for (int j = 0; j < Columns_; ++j)
        AX[i] += A_[i*Columns_+j]*X[j];
    }
    return AX;
  }
  
  /*!
    @brief  Matrix-Matrix multiply
    @param[in] B   matrix to post-multiply with
  */
  
  Matrix operator*(Matrix const& B) const {
    assert(Columns_ == B.rows());
    
    Matrix AB(Rows_, B.columns(), 0.0);
    
    for (int i = 0; i < Rows_; ++i)
      for (int j = 0; j < B.columns(); ++j)
        for (int k = 0; k < Columns_; ++k)
          AB[i][j] += A_[i*Columns_+k]*B[k][j];
    
    return AB;
  }

  /*!
    @brief  Matrix data
    @return   Flat vector with matrix entries
  */
  const std::vector<double>& data() const { return A_; }

  /*!
    @brief  solve a linear system A X = B with this matrix (A)
    @param[in] B  right-hand sides (multiple)
    @param[in] method what method to use for solution
    @param[in,out] error message, if any 
    @return the solution X

    method=="inverse" ==> use the  inverse operator
    method=="lapack-posv" ==> use lapack dposvx for symmetric positive definite A.
    method=="lapack-sysv" ==> use lapack dsysvx for symmetric A.
    method=="lapack-gesv" ==> use lapack dgesvx for general A.
    method=="lapack-sytr" ==> use lapack dsytrf & dsytrs for symmetric A.

    If \code error\endcode is not present or has value "ignore", no
    message will be returned.
    The value of \code error\endcode will be "ignore" on return.
    If \code error\endcode is present and has a value other than
    "ignore", the value "none" will be returned if no error was
    generated, or else will contain the appropriate error message.
  */
  Matrix solve(Matrix const& B,
               std::string method="inverse",
	       std::string &error=ignore_msg_ref);

  /*! 
    @brief report whether matrix is singular or not
    @return if 0: unknown, if 1: not singular, if 2: singular
    A value of 0 indicates the inverse() or solve() methods have not been run.
    A value of 1 indicates the methods have run and the result is definitely singular.
    A value of 2 indicates the methods have run and the result is definitely non-singular.
  */
  int is_singular(){
    return is_singular_;
  }

 private:
  int Rows_, Columns_;
  std::vector<double> A_;
  std::string method_;
  int is_singular_ = 0;
};  // class Matrix

// Add two matrices.
inline
const Matrix operator+(const Matrix& A, const Matrix& B) {
  assert((A.rows() == B.rows()) && (A.columns() == B.columns()));

  Matrix Sum(A);
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      Sum[i][j] += B[i][j];
  
  return Sum;
}

/// Subtract two materices.
inline
const Matrix operator-(const Matrix& A, const Matrix& B) {
  assert((A.rows() == B.rows()) && (A.columns() == B.columns()));

  Matrix Diff(A);
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      Diff[i][j] -= B[i][j];
  
  return Diff;
}

/*!
  @brief  Multiply a Matrix by a scalar
  @param[in] x The scaling factor 
*/
inline  
const Matrix operator*(const Matrix& A, const double& s) {
  Matrix As(A);
  
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      As[i][j] *= s;
  
  return As;
}

/*!
  @brief  Multiply a Matrix by a scalar
  @param[in] x The scaling factor
*/
inline  
const Matrix operator*(const double& s, const Matrix& A) {
  Matrix sA(A);
  
  for (int i = 0; i < A.rows(); ++i)
    for (int j = 0; j < A.columns(); ++j)
      sA[i][j] *= s;
  
  return sA;
}

// Multiply the first vector by the transpose of the second vector
template<int D>
inline
Matrix operator*(const Vector<D>& a, const Vector<D>& b) {
  Matrix prod(D, D);
  for (int i = 0; i < D; i++) 
    for (int j = 0; j < D; j++)
      prod[i][j] = a[i]*b[j];
  return prod;
}

/*!
  @brief  Computes the solution of A*x = b using the QR decomposition
  @tparam D Indicates the dimensionality of the Vector
  @param[in] A  The system matrix
  @param[in] b  The system right-hand side
  @param[out] x  The solution vector
*/
template<int D>
inline
void solve(const Matrix& A, const Vector<D>& b, Vector<D>& x) {
  int n = A.rows();
  assert(n == A.columns());
  assert(D == n);

  const std::vector<double>& a_ = A.data();
  double s, norm;
   
  Matrix R(A);
  Matrix Q(n, n, 0.0);
  Matrix U(n, n);

  std::vector<double> y(n);
  
  int i, j, k, l;
  for (i = 0; i < n; i++)
    Q[i][i] = 1.0;
  
  // Find A = QR using the reflection method
  for (i = 0; i < n - 1; i++) {
    y[i] = R[i][i];
    for (j = i + 1, s = 0.0; j < n; j++) {
      y[j] = R[j][i];
      s += pow2(R[j][i]);
    }
    norm = sqrt(s + pow2(y[i]));
    assert(std::fabs(norm > std::numeric_limits<double>::epsilon()));
    if (s == 0)
      continue;
    
    y[i] -= norm;
    norm = sqrt(s + pow2(y[i]));

    for (j = i; j < n; j++)
      y[j] /= norm;
    
    // Generate U matrix
    for (k = 0; k < n; k++)
      for (l = 0; l < n; l++)
        if (k < i || l < i)
          U[k][l] = (k == l) ? 1.0 : 0.0;
        else
          U[k][l] = (k == l) ? 1.0 - 2.0*y[k]*y[l] : -2.0*y[k]*y[l];
    
    // Modify Q matrix
    Q = Q*U;
    // Modify R matrix
    R = U*R;
  }
  // We need to confirm that the last element is nonzero
  assert(std::fabs(R[n - 1][n - 1]) > std::numeric_limits<double>::epsilon());

  Vector<D> QTb = Q.transpose()*b;
  
  // Now, use the back substitution to find x
  for (i = n - 1; i >= 0; i--) {
    for (k = i + 1; k < n; k++)
      QTb[i] -= R[i][k]*x[k];

    x[i] = QTb[i]/R[i][i];
  }
}

/*!
  @brief  Computes the solution of A*x = b for 1x1 algebraic system
  @param[in] A  The system matrix
  @param[in] b  The system right-hand side
  @param[out] x  The solution vector
*/
template<>
inline
void solve<1>(const Matrix& A, const Vector<1>& b, Vector<1>& x) {
  assert(std::fabs(A[0][0]) > std::numeric_limits<double>::epsilon());
  x[0] = b[0]/A[0][0];
}

/*!
  @brief  Computes the solution of A*x = b for 2x2 algebraic system
  @param[in] A  The system matrix
  @param[in] b  The system right-hand side
  @param[out] x  The solution vector
*/
template<>
inline
void solve<2>(const Matrix& A, const Vector<2>& b, Vector<2>& x) {
  double detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
  assert(std::fabs(detA) > std::numeric_limits<double>::epsilon());

  x[0] = (A[1][1]*b[0] - A[0][1]*b[1])/detA;
  x[1] = (A[0][0]*b[1] - A[1][0]*b[0])/detA;
}

}  // namespace Wonton

#endif  // WONTON_MATRIX_H_
