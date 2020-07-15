/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/



#ifndef WONTON_LS_FITS_H_
#define WONTON_LS_FITS_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Point.h"
#include "wonton/support/Matrix.h"
#include "wonton/support/svd.h"

namespace Wonton {


/**
 * @brief Build the matrix A and its transpose A^T used in the
 *        equation (A^T.A).X = (A^T.F) from the stencil.
 *
 * Since the stencil is the same for multiple field variables,
 * it should be called only once. It will also be used to compute
 * the right hand side of the equation. Here the first point
 * of the stencil is assumed to be the point where the gradient
 * must be computed.
 *
 * @tparam dim: the spatial dimension of the problem.
 * @param stencil: the stencil point coordinates.
 * @return the matrices A^T and A.
 */
template<int dim>
std::pair<Matrix, Matrix> build_transposed_matrix(std::vector<Point<dim>> const& stencil) {

  // the first point is the reference point where we are trying
  // to compute the gradient; so the number of rows is size - 1.
  int const num = stencil.size() - 1;

  // each row of A contains the components of the vector from
  // coords[0] to the candidate point coords[i] being used
  // in the least squares approximation (x_i - x_0).
  Matrix A(num, dim);
  for (int i = 0; i < num; ++i) {
    for (int j = 0; j < dim; ++j) {
      A[i][j] = stencil[i+1][j] - stencil[0][j];
    }
  }

  return std::make_pair(A.transpose(), A);
}

/**
 * @brief Build the right hand side (A^T.F) of the least square
 *        equation (A^T.A).X = (A^T.F) from field values.
 *
 * It should be called for each field variable using the cached
 * transposed matrix A^T used in the least square equation.
 *
 * @tparam dim: the spatial dimension of the problem.
 * @param AT: the cached transposed matrix A^T.
 * @param values: the field value at each stencil point.
 * @return the vector (A^T.F)
 */
template<int dim>
Vector<dim> build_right_hand_side(Matrix const& AT, std::vector<double> const& values) {

  // the first value is that of the reference point where we are trying
  // to compute the gradient; so the number of components is size - 1.
  int const num = values.size() - 1;

  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)
  Vector<num> F;
  for (int i = 0; i < num; ++i) {
    F[i] = values[i+1] - values[0];
  }

  return AT * F;
}

/**
 * @brief Compute the gradient using the cached inverse of the
 *        (A^T.A) matrix and the given right hand side (A^T.F).
 *
 * @tparam dim: the spatial dimension of the problem.
 * @tparam CoordSys: the current coordinate system.
 * @param ATA_inv: the cached inverse of the least square matrix.
 * @param ATF: the right hand side of the least square equation.
 * @param reference: the coordinate of the reference point.
 * @return the gradient at the current point.
 */
template<int dim, typename CoordSys = CartesianCoordinates>
Vector<dim> ls_gradient(Matrix const& ATA_inv,
                        Vector<dim> const& ATF,
                        Point<dim> const& reference) {

  CoordSys::template verify_coordinate_system<dim>();

  // compute gradient
  Vector<dim> gradient = ATA_inv * ATF;

  // correct it for curvilinear coordinates
  CoordSys::modify_gradient(gradient, reference);

  return gradient;
}

/*!
  @brief Compute least squares gradient from set of values
  @param[in] coords Vector of coordinates at which values are given
  @param[in] vals   Vector of values at said coordinates

  Compute a least squares gradient from a set of values. The first
  point is assumed to be the point where the gradient must be computed
  and the first value is assumed to the value at this reference point

  This operator does not know anything about a mesh.
  It solves the algebraic equation using an optimized LAPACK kernel
  if available, or the inverse method if not. It is meant to be used
  only when remapping a single field variable.

*/

template<int D, typename CoordSys = CartesianCoordinates>
Vector<D> ls_gradient(std::vector<Point<D>> const & coords,
                      std::vector<double> const & vals) {

  CoordSys::template verify_coordinate_system<D>();

  // construct the least square equation components
  auto const M = build_transposed_matrix(coords);
  Matrix ATA = M.first * M.second;
  Vector<D> ATF = build_right_hand_side(M.first, vals);

  // solve it
#ifdef WONTON_HAS_LAPACKE
  // use lapack solver for symmetric positive definite matrices
  Vector<D> gradient = ATA.solve(ATF, "lapack-posv");
#else
  // use the inverse method
  Vector<D> gradient = ATA.solve(ATF, "inverse");
#endif

  // correct it for curvilinear coordinates
  CoordSys::modify_gradient(gradient, coords[0]);

  return gradient;
}


/*!
  @brief Compute least squares quadfit from set of values
  @param[in] coords Vector of coordinates at which values are given
  @param[in] vals   Vector of values at said coordinates

  Compute a least squares quadfit from a set of values. The first
  point is assumed to be the point where the quadfit must be computed
  and the first value is assumed to the value at this reference point

  This operator does not know anything about a mesh.

*/

template<int D>
// int N = D*(D+3)/2;
Vector<D*(D+3)/2> ls_quadfit(std::vector<Point<D>> const & coords,
                             std::vector<double> const & vals,
                             bool const boundary_element) {

  double val0 = vals[0];

  // There are nvals but the first is the reference point where we
  // are trying to compute the quadfit; so the matrix sizes etc
  // will only be nvals-1

  int nvals = vals.size();

  if (boundary_element) {
    // Not enough values to do a quadratic fit - drop down to linear fit

    Vector<D> grad = ls_gradient(coords, vals);
    Vector<D*(D+3)/2> result;
    for (int i = 0; i < D*(D+3)/2; i++) result[i] = 0.0;
    // Fill in the begining of the vector with gradient terms
    for (int i = 0; i < D; i++) result[i] = grad[i];
    return result;
  }


  // Each row of A contains the components of the vector from
  // coord0 to the candidate point being used in the Least Squares
  // approximation (X_i-X_0), up to quadratic order. Index j0
  // labels the columns up to D (e.g., gradient). Index j1
  // labels the D(D+1)/2 extra columns for the quadratic terms.
  // Index i labels rows for data points away from the center of
  // the point cloud (i.e., array[0]) with nvals points, contained
  // in the array coords(nvals,D).

  Matrix A(nvals-1, D*(D+3)/2); // don't include 0th point
  for (int i = 0; i < nvals-1; i++) {
    int j1 = D; // index of colunns with quadrtic terms
    for (int j0 = 0; j0 < D; ++j0) {
      A[i][j0] = coords[i+1][j0]-coords[0][j0];
      // Add columns with the remaining quadratic delta terms
      // for each i, starting at D+j, reusing linear delta terms
      for (int k=0; k <= j0; k++) {
        A[i][j1] = A[i][k]*A[i][j0];
        j1 += 1;
      }
    }
  }

  int const ma = D*(D+3)/2;
  std::vector<std::vector<double> > u(nvals-1, std::vector<double>(ma));
  std::vector<std::vector<double> > v(ma, std::vector<double>(ma));
  std::vector<double> w(ma);
  for (int i = 0; i < nvals-1; i++) {
    for (int j = 0; j < ma; j++) {
      u[i][j] = A[i][j];
    }
  }

  svd(u,w,v);

  // "edit" the singular values (eigenvalues):
  // (1) find wmax = largest value of w
  // (2) select only w's between TOL*wmax <= w <= wmax
  // -->Any other w's contribute 0 to the solution
  double scale = 1e-5;
  double wmax = 0.0;
  for (int j = 0; j < ma; j++) {
    if (w[j] > wmax) wmax=w[j];
  }
  double thresh = scale*wmax;
  for (int j = 0; j < ma; j++) {
    if (w[j] < thresh) w[j]=0.0;
  }

  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)

  std::vector<double> F(nvals-1);
  for (int i = 0; i < nvals-1; ++i) {
    F[i] = vals[i+1]-val0;
  }

  // solve the problem for coefficients, in "result".
  std::vector<double> result(ma);
  svd_solve(u,w,v,F,result);

  return result;
}

}  // namespace Wonton

#endif
