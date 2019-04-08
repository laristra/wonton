/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_DIRECT_PRODUCT_MESH_H_
#define WONTON_DIRECT_PRODUCT_MESH_H_

#include <cassert>
#include <vector>

#include "wonton/support/wonton.h"

/*!
  @file direct_product_mesh.h
  @brief Definition of the Direct_Product_Mesh class.

  A Direct_Product_Mesh is a basic, serial, N-dimensional, axis-aligned,
  logically-rectangular mesh.  It is called a direct product mesh because it is
  the direct product of independent discretizations along each axis.
 */

namespace Wonton {

/*!
  @class Direct_Product_Mesh "direct_product_mesh.h"
  @brief A basic, axis-aligned, logically-rectangular mesh.

  A Direct_Product_Mesh is a basic, serial, N-dimensional, axis-aliged,
  logically-rectangular mesh.  It is called a direct product mesh because it is
  the direct product of independent discretizations along each axis.

  The discretizations need not be uniform -- that is, the cell sizes can vary
  across the mesh.  However, they are static.  Once set, they will not change.

  Storage of the mesh is based on cell edge coordinates, aka points.

  The Direct_Product_Mesh is designed to implement only the necessary
  functionality to test certain components in Wonton and Portage.  As the scope
  of the tests expands, the scope of functionality of the Direct_Product_Mesh
  may also expand.
 */
template<int D>
class Direct_Product_Mesh {

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  Direct_Product_Mesh() = delete;

  /*!
    @brief Constructor to create a Direct_Product_Mesh.
    @param[in] edges The cell edge coordinates along each axis.

    Specify the edge coordinates that delineate the cells of the mesh along
    each axis.
  */
  Direct_Product_Mesh(const std::vector<double> edges[D]);

  //! Assignment operator (disabled).
  Direct_Product_Mesh & operator=(const Direct_Product_Mesh<D> &) = delete;

  //! Destructor
  ~Direct_Product_Mesh();


  // ==========================================================================
  // Accessors

  /*!
    @brief Get the dimensionality of the mesh.
  */
  int space_dimension() const;

  /*!
    @brief Get the number of points (edge coordinates) along a specified axis.
    @param[in] dim The axis to query.
  */
  int axis_num_points(const int dim) const;

  /*!
    @brief Get the specified point (edge coordinate).
    @param[in] dim The axis to query.
    @param[in] dim The point to query.
  */
  double axis_point_coordinate(const int dim, const int pointid) const;


 private:

  // ==========================================================================
  // Class data

  //! Cell edge coordinates along the three axes
  std::vector<double> edges_[D];

};  // class Direct_Product_Mesh


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor
template<int D>
Direct_Product_Mesh<D>::Direct_Product_Mesh(
    const std::vector<double> edges[D]) {
  for (int d = 0; d < D; ++d) {
    edges_[d] = edges[d];
  }
}

// ____________________________________________________________________________
// Destructor
template<int D>
Direct_Product_Mesh<D>::~Direct_Product_Mesh() {
  for (int d = 0; d < D; ++d) {
    edges_[d].clear();
  }
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get the dimensionality of the mesh.
template<int D>
int Direct_Product_Mesh<D>::space_dimension() const {
  return(D);
}

// ____________________________________________________________________________
// Get the number of points (edge coordinates).
template<int D>
int Direct_Product_Mesh<D>::axis_num_points(const int dim) const {
  assert(dim >= 0);
  assert(dim < D);
  return edges_[dim].size();
}

// ____________________________________________________________________________
// Get the specified point (edge coordinate).
template<int D>
double Direct_Product_Mesh<D>::axis_point_coordinate(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);
  assert(pointid >= 0);
  assert(pointid < edges_[dim].size());
  return edges_[dim][pointid];
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
