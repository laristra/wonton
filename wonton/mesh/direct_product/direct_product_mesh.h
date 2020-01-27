/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_DIRECT_PRODUCT_MESH_H_
#define WONTON_DIRECT_PRODUCT_MESH_H_

#include <array>
#include <cassert>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/CoordinateSystem.h"

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

  Storage of the mesh is based on cell boundary coordinates, aka axis points.

  The Direct_Product_Mesh is designed to implement only the necessary
  functionality to test certain components in Wonton and Portage.  As the scope
  of the tests expands, the scope of functionality of the Direct_Product_Mesh
  may also expand.
 */
template<int D, class CoordSys=Wonton::DefaultCoordSys>
class Direct_Product_Mesh {

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  Direct_Product_Mesh() = delete;

  /*!
    @brief Constructor to create a Direct_Product_Mesh.
    @param[in] axis_points The cell boundary coordinates along each axis.

    Specify the axis point coordinates that delineate the cells of the mesh
    along each axis.
  */
  Direct_Product_Mesh(const std::array<std::vector<double>,D> &axis_points);

  //! Assignment operator (disabled).
  Direct_Product_Mesh & operator=(
      const Direct_Product_Mesh<D,CoordSys> &) = delete;

  //! Destructor
  ~Direct_Product_Mesh();


  // ==========================================================================
  // Accessors

  /*!
    @brief Get the dimensionality of the mesh.
  */
  int space_dimension() const;

  /*!
    @brief Get the number of axis points (cell boundary coordinates) along a
    specified axis.  @param[in] dim The axis to query.
  */
  int num_axis_points(const int dim) const;

  /*!
    @brief Get the specified axis point (cell boundary coordinate).
    @param[in] dim The axis to query.
    @param[in] dim The point to query.
  */
  double get_axis_point(const int dim, const int pointid) const;


 private:

  // ==========================================================================
  // Class data

  //! Axis points (cell boundary coordinates) along the three axes
  std::array<std::vector<double>,D> axis_points_;

};  // class Direct_Product_Mesh


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor
template<int D, class CoordSys>
Direct_Product_Mesh<D,CoordSys>::Direct_Product_Mesh(
    const std::array<std::vector<double>,D> &axis_points) {
  for (int d = 0; d < D; ++d) {
    axis_points_[d] = axis_points[d];
  }
}

// ____________________________________________________________________________
// Destructor
template<int D, class CoordSys>
Direct_Product_Mesh<D,CoordSys>::~Direct_Product_Mesh() {
  for (int d = 0; d < D; ++d) {
    axis_points_[d].clear();
  }
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get the dimensionality of the mesh.
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::space_dimension() const {
  return(D);
}

// ____________________________________________________________________________
// Get the number of axis points (cell boundary coordinates).
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::num_axis_points(const int dim) const {
  assert(dim >= 0);
  assert(dim < D);
  return axis_points_[dim].size();
}

// ____________________________________________________________________________
// Get the specified axis point (cell boundary coordinate).
template<int D, class CoordSys>
double Direct_Product_Mesh<D,CoordSys>::get_axis_point(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);
  assert(pointid >= 0);
  assert(pointid < axis_points_[dim].size());
  return axis_points_[dim][pointid];
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
