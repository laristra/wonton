/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
#define WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_

#include "wonton/mesh/direct_product/direct_product_mesh.h"

#include <vector>
#include <algorithm>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

/*!
  @file direct_product_mesh_wrapper.h
  @brief Wrapper for a Direct_Product_mesh
 */

namespace Wonton {
  /*!
    @class Direct_Product_Mesh_Wrapper direct_product_mesh_wrapper.h
    @brief A thin wrapper that implements mesh methods for Direct_Product_Mesh

    The methods implemented are those required by select parts of the Wonton
    and Portage code.  This will expand as the list of components that this
    wrapper supports expands.
   */
class Direct_Product_Mesh_Wrapper {

 public:

  // ==========================================================================
  // Constructors and destructors

  /*!
    @brief Constructor for the mesh wrapper.
    @param[in] mesh The Direct_Product_Mesh we wish to wrap.
   */
  explicit Direct_Product_Mesh_Wrapper(Direct_Product_Mesh const & mesh);

  //! Copy constructor (disabled).
  Direct_Product_Mesh_Wrapper(Direct_Product_Mesh_Wrapper const &) = delete;

  //! Assignment operator (disabled).
  Direct_Product_Mesh_Wrapper & operator=(
      Direct_Product_Mesh_Wrapper const &) = delete;

  //! Destructor
  ~Direct_Product_Mesh_Wrapper();

  // ==========================================================================
  // Accessors

  // TODO: doxygen
  template<long D>
  void get_global_bounds(Point<D> *plo, Point<D> *phi) const;

  // TODO: doxygen
  counting_iterator axis_point_begin(const int dim) const;

  // TODO: doxygen
  counting iterator axis_point_end(const int dim) const;

  // TODO: doxygen
  double axis_point_coordinate(const int dim, const int pointid) const;

  // TODO: doxygen
  int axis_num_cells(const int dim) const;


 private:

  // ==========================================================================
  // Class data

  //! The mesh to wrap.
  Direct_Product_Mesh const & mesh_;
};  // class Direct_Product_Mesh_Wrapper


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// constructor
explicit Direct_Product_Mesh_Wrapper::Direct_Product_Mesh_Wrapper(
    Direct_Product_Mesh const & mesh) :
    mesh_(mesh) {
}

// ____________________________________________________________________________
// destructor
Direct_Product_Mesh_Wrapper::~Direct_Product_Mesh_Wrapper() {
}

// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get global mesh bounds
template<long D>
void get_global_bounds(Point<D> *plo, Point<D> *phi) const {
  assert(D == mesh_.space_dimension())
  (*plo)[0] = edges_i_.front();
  (*phi)[0] = edges_i_.back();
  if (D >= 2) {
    (*plo)[1] = edges_j_.front();
    (*phi)[1] = edges_j_.back();
    if (D >= 3) {
      (*plo)[2] = edges_k_.front();
      (*phi)[2] = edges_k_.back();
    }
  }
}

// ____________________________________________________________________________
// Get iterator for axis edge coordinates (beginning of array)
counting_iterator axis_point_begin(const int dim) const {
  int start_index = 0;
  return make_counting_iterator(start_index);
}


// ____________________________________________________________________________
// Get iterator for axis edge coordinates (end of array)
counting iterator axis_point_end(const int dim) const {
  int start_index = 0;
  return make_counting_iterator(start_index + mesh.axis_num_points(dim));
}


// ____________________________________________________________________________
// Get edge coordinate value
double axis_point_coordinate(const int dim, const int pointid) const {
  return mesh.axis_point_coordinate(dim, pointid);
}


// ____________________________________________________________________________
// Get number of cells along axis
int axis_num_cells(const int dim) const {
  return mesh.axis_num_points(dim) - 1;
}


}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
