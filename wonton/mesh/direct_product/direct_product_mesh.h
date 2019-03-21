/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_DIRECT_PRODUCT_MESH_H_
#define WONTON_DIRECT_PRODUCT_MESH_H_

#include <vector>
#include <memory>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

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

  The Direct_Product_Mesh is designed to implement only the necessary
  functionality to test certain components in Wonton and Portage.  As the scope
  of the tests expands, the scope of functionality of the Direct_Product_mesh
  may also expand.
 */
class Simple_Mesh {

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

  //! Dimensionality of the mesh
  int dimensionality_ = -1;

  //! Cell edge coordinates along the three axes
  std::vector<double> edges_i_, edges_j_, edges_k_;

};  // class Direct_Product_Mesh

// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// 1D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_mesh(
    const std::vector<double> edges) :
    edges_i_(edges) {
  dimensionality_ = 1;
  set_default_coordinates();
}

// ____________________________________________________________________________
// 2D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_mesh(
    const std::vector<double> edges_i,
    const std::vector<double> edges_j) :
    edges_i_(edges_i),
    edges_j_(edges_j) {
  // TODO: Note that this may have to be modified.  For example, in 2D you can
  // have cylindrical coordinates (s,z) or polar coordinates (s, phi), and
  // fitting those both into the standard 3D cylindrical (s, phi, z) requires
  // that 2D cylindrical will assign edges_j to edges_k_ not edges_j_.
  dimensionality_ = 2;
  set_default_coordinates();
}

// ____________________________________________________________________________
// 3D constructor (see class definition for declaration).
Direct_Product_Mesh::Direct_Product_mesh(
    const std::vector<double> edges_i,
    const std::vector<double> edges_j,
    const std::vector<double> edges_k) :
    edges_i_(edges_i),
    edges_j_(edges_j),
    edges_j_(edges_k) {
  dimensionality_ = 3;
  set_default_coordinates();
}

// ____________________________________________________________________________
// destructor (see class definition for declaration).
Direct_Product_Mesh::~Direct_Product_mesh() {
  dimensionality_ = -1;
  edges_i_.clear();
  edges_j_.clear();
  edges_k_.clear();
}


// ============================================================================
// Private support methods
void set_default_coordinates() {
  switch(dimensionality_) {
    case 1 :
      // Cartesian coordinates (currently no others available)
      edges_j_ = { 0.0, 1.0 };
      edges_k_ = { 0.0, 1.0 };
      break;
    case 2 :
      // Cartesian coordinates (currently no others available)
      edges_k_ = { 0.0, 1.0 };
      break;
    case 3 :
      // All edge coordinate arrays are specified, so do nothing.  This clause
      // exists to avoid the default clause.
      break;
    default :
      // TODO: Generate an error.
  }
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
