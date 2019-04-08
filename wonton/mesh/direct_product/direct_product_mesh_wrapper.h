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
#include "wonton/support/CellID.h"
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

  //! Get iterator for axis edge coordinates (beginning of array).
  counting_iterator axis_point_begin(const int dim) const;

  //! Get iterator for axis edge coordinates (end of array).
  counting_iterator axis_point_end(const int dim) const;

  //! Get edge coordinate value.
  double axis_point_coordinate(const int dim, const int pointid) const;

  //! Get number of cells along axis.
  int axis_num_cells(const int dim) const;

  //! Get number of cells in entire mesh.
  int total_num_cells() const;

  //! Get lower and upper corners of cell bounding box
  void cell_get_bounds(const CellID id, Point<D> *plo, Point<D> *phi) const;

  // ==========================================================================
  // Index/ID conversions

  //! Convert from indices to cell ID
  CellID indices_to_cellid(const std::array<int,D>& indices) const;

  //! Convert from ID to indices
  std::array<int,D> cellid_to_indices(const CellID id) const;


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
    (*plo)[d] = mesh_.axis_point_coordinate(d,0);
    (*phi)[d] = mesh_.axis_point_coordinate(d,mesh_.axis_num_points(d)-1);
  }
}

// ____________________________________________________________________________
// Get iterator for axis edge coordinates (beginning of array).
template<int D>
counting_iterator Direct_Product_Mesh_Wrapper<D>::axis_point_begin(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  int start_index = 0;
  return make_counting_iterator(start_index);
}

// ____________________________________________________________________________
// Get iterator for axis edge coordinates (end of array).
template<int D>
counting_iterator Direct_Product_Mesh_Wrapper<D>::axis_point_end(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  int start_index = 0;
  return make_counting_iterator(start_index + mesh_.axis_num_points(dim));
}

// ____________________________________________________________________________
// Get edge coordinate value.
template<int D>
double Direct_Product_Mesh_Wrapper<D>::axis_point_coordinate(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.axis_point_coordinate(dim, pointid);
}

// ____________________________________________________________________________
// Get number of cells along axis.
template<int D>
int Direct_Product_Mesh_Wrapper<D>::axis_num_cells(const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.axis_num_points(dim) - 1;
}

// ____________________________________________________________________________
// Get number of cells in entire mesh.
template<int D>
int Direct_Product_Mesh_Wrapper<D>::total_num_cells() const {
  int count = 1;
  for (int dim = 0; dim < mesh_.space_dimension(); ++dim) {
    count *= mesh_.axis_num_points(dim) - 1;
  }
  return count;
}

// ____________________________________________________________________________
// Get lower and upper corners of cell bounding box
template<int D>
void Direct_Product_Mesh_Wrapper<D>::cell_get_bounds(
    const CellID id, Point<D> *plo, Point<D> *phi) const {
  std::array<int,D> indices = cellid_to_indices(id);
  // Cell edges (points) are zero-indexed, cells are zero-indexed.  Thus cell 0
  // is bounded by edges 0 and 1, or more generally, cell N is bounded by edges
  // N and N+1.
  for (int d = 0; d < D; ++d) {
    (*plo)[d] = mesh_.axis_point_coordinate(d, indices[d]);
    (*phi)[d] = mesh_.axis_point_coordinate(d, indices[d]+1);
  }
}

// ============================================================================
// Index/ID conversions

// ____________________________________________________________________________
// Convert from indices to ID
template<int D>
CellID Direct_Product_Mesh_Wrapper<D>::indices_to_cellid(
    const std::array<int,D>& indices) const {
  assert(D == mesh_.space_dimension());
  CellID id = 0;
  switch(D) {
    case 3 :
    {
      CellID k = (CellID) indices[2];
      CellID jmax = (CellID) axis_num_cells(1);
      id += k;
      id *= jmax;
    }
    case 2 :
    {
      CellID j = (CellID) indices[1];
      CellID imax = (CellID) axis_num_cells(0);
      id += j;
      id *= imax;
    }
    case 1 : 
    {
      CellID i = (CellID) indices[0];
      id += i;
    }
  }
  return id;
}

// ____________________________________________________________________________
// Convert from ID to indices
template<int D>
std::array<int,D> Direct_Product_Mesh_Wrapper<D>::cellid_to_indices(
    const CellID id) const {
  assert(D == mesh_.space_dimension());
  std::array<int,D> indices;
  CellID index, denom;
  CellID residual = id;
  switch(D) {
    case 3 :
      denom = (CellID) (axis_num_cells(1) * axis_num_cells(0));
      index = residual / denom;
      residual -= index * denom;
      indices[2] = (int) index;
    case 2 :
      denom = (CellID) axis_num_cells(0);
      index = residual / denom;
      residual -= index * denom;
      indices[1] = (int) index;
    case 1 :
      indices[0] = (int) residual;
  }
  return std::move(indices);
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
