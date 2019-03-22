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

  A Direct_Product_Mesh is a basic, serial, 1/2/3D, axis-aligned,
  logically-rectangular mesh.  It is called a direct product mesh because it is
  the direct product of independent discretizations along each axis.
 */

namespace Wonton {

/*!
  @class Direct_Product_Mesh "direct_product_mesh.h"
  @brief A basic, axis-aligned, logically-rectangular mesh.

  A Direct_Product_Mesh is a basic, serial, 1D/2D/3D, axis-aliged,
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
class Direct_Product_Mesh {

 public:

  // ==========================================================================
  // Constructors and destructors

  /*!
    @brief Constructor to create a 1D Direct_Product_Mesh.
    @param[in] edges The cell edge coordinates.

    Specify the edge coordinates that delineate the cells of the mesh.
  */
  Direct_Product_Mesh(std::vector<double> edges);

  /*!
    @brief Constructor to create a 2D Direct_Product_Mesh.
    @param[in] edges_i The cell edge coordinates along the first axis.
    @param[in] edges_j The cell edge coordinates along the second axis.

    Specify the edge coordinates that delineate the cells of the mesh along
    each axis.
  */
  Direct_Product_Mesh(std::vector<double> edges_i,
                      std::vector<double> edges_j);

  /*!
    @brief Constructor to create a 3D Direct_Product_Mesh.
    @param[in] edges_i The cell edge coordinates along the first axis.
    @param[in] edges_j The cell edge coordinates along the second axis.
    @param[in] edges_k The cell edge coordinates along the third axis.

    Specify the edge coordinates that delineate the cells of the mesh along
    each axis.
  */
  Direct_Product_Mesh(std::vector<double> edges_i,
                      std::vector<double> edges_j,
                      std::vector<double> edges_k);

  //! Assignment operator (disabled).
  Direct_Product_Mesh & operator=(const Direct_Product_Mesh &) = delete;

  //! Destructor
  ~Direct_Product_Mesh();


  // ==========================================================================
  // Accessors

  /*!
    @brief Get the dimensionality of the mesh.
  */
  int space_dimension() const;

  /*!
    @brief Get the number of points (edge coordinates).
  */
  int axis_num_points(const int dim) const;

  /*!
    @brief Get the specified point (edge coordinate).
  */
  double axis_point_coordinate(const int dim, const int pointid) const;


 private:

  // ==========================================================================
  // Private support methods

  /*!
    @brief Constructs the default coordinates

    For unused axes (e.g., the third axis for a 2D mesh), create the cell edge
    coordinate arrays.

    This method is primarily intended for extensions to curvilinear
    coordinates.  In Cartesian coordinates, the axes are infinite, so we take a
    unit-length slice along that axis.  In curvilinear coordinates, this will
    depend on the coordinate system.
   */
  void set_default_coordinates();


  // ==========================================================================
  // Class data

  //! Maximum dimensionality allowed
  static const int MAXDIM_ = 3;

  //! Dimensionality of the mesh
  int dimensionality_ = -1;

  //! Cell edge coordinates along the three axes
  std::vector<double> edges_[MAXDIM_];

};  // class Direct_Product_Mesh


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// 1D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_Mesh(const std::vector<double> edges) {
  dimensionality_ = 1;
  edges_[0] = edges;
  set_default_coordinates();
}

// ____________________________________________________________________________
// 2D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_Mesh(
    const std::vector<double> edges_i,
    const std::vector<double> edges_j) {
  dimensionality_ = 2;
  edges_[0] = edges_i;
  edges_[1] = edges_j;
  set_default_coordinates();
}

// ____________________________________________________________________________
// 3D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_Mesh(
    const std::vector<double> edges_i,
    const std::vector<double> edges_j,
    const std::vector<double> edges_k) {
  dimensionality_ = 3;
  edges_[0] = edges_i;
  edges_[1] = edges_j;
  edges_[2] = edges_k;
  set_default_coordinates();
}

// ____________________________________________________________________________
// destructor (see class definition for declaration).
Direct_Product_Mesh::~Direct_Product_Mesh() {
  dimensionality_ = -1;
  for (int d = 0; d < MAXDIM_; ++d) {
    edges_[d].clear();
  }
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get the dimensionality of the mesh.
int Direct_Product_Mesh::space_dimension() const {
  return dimensionality_;
}

// ____________________________________________________________________________
// Get the number of points (edge coordinates).
int Direct_Product_Mesh::axis_num_points(const int dim) const {
  return edges_[dim].size();
}

// ____________________________________________________________________________
// Get the specified point (edge coordinate).
double Direct_Product_Mesh::axis_point_coordinate(
    const int dim, const int pointid) const {
  return edges_[dim][pointid];
}


// ============================================================================
// Private support methods
void Direct_Product_Mesh::set_default_coordinates() {
  assert(dimensionality_ >= 1);
  assert(dimensionality_ <= 3);
  switch(dimensionality_) {
    case 1 :
      // Cartesian coordinates (currently no others available)
      edges_[1] = { 0.0, 1.0 };
      edges_[2] = { 0.0, 1.0 };
      break;
    case 2 :
      // Cartesian coordinates (currently no others available)
      edges_[2] = { 0.0, 1.0 };
      break;
    case 3 :
      // All edge coordinate arrays are specified, so do nothing.  This clause
      // exists to avoid the default clause.
      break;
  }
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
