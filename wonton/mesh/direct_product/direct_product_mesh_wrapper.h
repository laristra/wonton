/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
#define WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_

#include "wonton/mesh/direct_product/direct_product_mesh.h"

#include <algorithm>
#include <array>
#include <vector>

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

  The methods implemented are those required by select parts of the Wonton and
  Portage code.  This will expand as the list of components that this wrapper
  supports expands.

  This class is not meant to be used to be utilized with downstream classes
  that assume an unstructured mesh (e.g. SearchKDTree or IntersectR2D etc.).
*/
template<int D>
class Direct_Product_Mesh_Wrapper {

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  Direct_Product_Mesh_Wrapper() = delete;

  /*!
    @brief Constructor for the mesh wrapper.
    @param[in] mesh The Direct_Product_Mesh we wish to wrap.
  */
  explicit Direct_Product_Mesh_Wrapper(Direct_Product_Mesh<D> const & mesh);

  //! Copy constructor (disabled).
  Direct_Product_Mesh_Wrapper(Direct_Product_Mesh_Wrapper<D> const &) = delete;

  //! Assignment operator (disabled).
  Direct_Product_Mesh_Wrapper & operator=(
      Direct_Product_Mesh_Wrapper<D> const &) = delete;

  //! Destructor
  ~Direct_Product_Mesh_Wrapper();

  // ==========================================================================
  // Accessors

  //! Get dimensionality of the mesh.
  int space_dimension () const;

  /*!
    @brief Get global mesh bounds.

    Because a Direct_Product_Mesh is required to be an axis-aligned box, the
    global mesh bounds are specified as a Point at the lower corner and a Point
    at the upper corner (lower and upper along all axes).

    The points must have the same size as the mesh dimensionality.
  */
  void get_global_bounds(Point<D> *plo, Point<D> *phi) const;

  //! Get iterator for axis point (beginning of array).
  counting_iterator axis_point_begin(const int dim) const;

  //! Get iterator for axis point (end of array).
  counting_iterator axis_point_end(const int dim) const;

  //! Get axis point value.
  double get_axis_point(const int dim, const int pointid) const;

  //! Get number of cells along axis.
  int axis_num_cells(const int dim) const;

  //! Get number of cells in entire mesh.
  int total_num_cells() const;

  //! Get lower and upper corners of cell bounding box
  void cell_get_bounds(const int id, Point<D> *plo, Point<D> *phi) const;

  // ==========================================================================
  // Index/ID conversions

  //! Convert from indices to cell ID
  int indices_to_cellid(const std::array<int,D>& indices) const;

  //! Convert from ID to indices
  std::array<int,D> cellid_to_indices(const int id) const;


 private:

  // ==========================================================================
  // Class data

  //! The mesh to wrap.
  Direct_Product_Mesh<D> const & mesh_;

};  // class Direct_Product_Mesh_Wrapper


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// constructor
template<int D>
Direct_Product_Mesh_Wrapper<D>::Direct_Product_Mesh_Wrapper(
    Direct_Product_Mesh<D> const & mesh) :
    mesh_(mesh) {
}

// ____________________________________________________________________________
// destructor
template<int D>
Direct_Product_Mesh_Wrapper<D>::~Direct_Product_Mesh_Wrapper() {
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get dimensionality of the mesh.
template<int D>
int Direct_Product_Mesh_Wrapper<D>::space_dimension() const {
  return mesh_.space_dimension();
}

// ____________________________________________________________________________
// Get global mesh bounds.
template<int D>
void Direct_Product_Mesh_Wrapper<D>::get_global_bounds(
    Point<D> *plo, Point<D> *phi) const {
  assert(D == mesh_.space_dimension());
  for (int d = 0; d < D; ++d) {
    (*plo)[d] = mesh_.get_axis_point(d,0);
    (*phi)[d] = mesh_.get_axis_point(d,mesh_.num_axis_points(d)-1);
  }
}

// ____________________________________________________________________________
// Get iterator for axis points (beginning of array).
template<int D>
counting_iterator Direct_Product_Mesh_Wrapper<D>::axis_point_begin(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  int start_index = 0;
  return make_counting_iterator(start_index);
}

// ____________________________________________________________________________
// Get iterator for axis points (end of array).
template<int D>
counting_iterator Direct_Product_Mesh_Wrapper<D>::axis_point_end(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  int start_index = 0;
  return make_counting_iterator(start_index + mesh_.num_axis_points(dim));
}

// ____________________________________________________________________________
// Get axis point value.
template<int D>
double Direct_Product_Mesh_Wrapper<D>::get_axis_point(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.get_axis_point(dim, pointid);
}

// ____________________________________________________________________________
// Get number of cells along axis.
template<int D>
int Direct_Product_Mesh_Wrapper<D>::axis_num_cells(const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.num_axis_points(dim) - 1;
}

// ____________________________________________________________________________
// Get number of cells in entire mesh.
template<int D>
int Direct_Product_Mesh_Wrapper<D>::total_num_cells() const {
  int count = 1;
  for (int dim = 0; dim < mesh_.space_dimension(); ++dim) {
    count *= mesh_.num_axis_points(dim) - 1;
  }
  return count;
}

// ____________________________________________________________________________
// Get lower and upper corners of cell bounding box
template<int D>
void Direct_Product_Mesh_Wrapper<D>::cell_get_bounds(
    const int id, Point<D> *plo, Point<D> *phi) const {
  std::array<int,D> indices = cellid_to_indices(id);
  // Cell axis points are zero-indexed, cells are zero-indexed.  Thus cell 0 is
  // bounded by axis points 0 and 1, or more generally, cell N is bounded by
  // axis points N and N+1.
  for (int d = 0; d < D; ++d) {
    (*plo)[d] = mesh_.get_axis_point(d, indices[d]);
    (*phi)[d] = mesh_.get_axis_point(d, indices[d]+1);
  }
}

// ============================================================================
// Index/ID conversions

// ____________________________________________________________________________
// Convert from indices to ID
template<int D>
int Direct_Product_Mesh_Wrapper<D>::indices_to_cellid(
    const std::array<int,D>& indices) const {
  assert(D == mesh_.space_dimension());
  int id = 0;
  // Loop over all but the last dimension
  for (int d = D-1; d > 0; --d) {
    int idx = indices[d];
    int mult = axis_num_cells(d-1);
    id += idx;
    id *= mult;
  }
  // Handle the last dimension
  id += indices[0];
  // Return
  return id;
}

// ____________________________________________________________________________
// Convert from ID to indices
template<int D>
std::array<int,D> Direct_Product_Mesh_Wrapper<D>::cellid_to_indices(
    const int id) const {
  assert(D == mesh_.space_dimension());
  std::array<int,D> indices;
  int residual = id;
  // Construct the denominators
  std::array<int,D> denom;
  denom[0] = 1;
  for (int d = 1; d < D; ++d) {
    denom[d] = denom[d-1] * axis_num_cells(d-1);
  }
  // Loop over all but the last dimension
  for (int d = D-1; d > 0; --d) {
    int index = residual / denom[d];
    residual -= index * denom[d];
    indices[d] = (int) index;
  }
  // Handle the last dimension
  indices[0] = (int) residual;
  // Return
  return std::move(indices);
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
