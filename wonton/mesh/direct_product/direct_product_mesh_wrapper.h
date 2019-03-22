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
#include "wonton/support/CellID.h"
#include "wonton/support/IntPoint.h"
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

  //! Get dimensionality of the mesh.
  int space_dimension () const;

  /*!
    @brief Get global mesh bounds.

    Because a Direct_Product_Mesh is required to be an axis-aligned box, the
    global mesh bounds are specified as a Point at the lower corner and a Point
    at the upper corner (lower and upper along all axes).

    The points must have the same size as the mesh dimensionality.
  */
  template<long D>
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
  int total_num_cells(const int dim) const;

  // ==========================================================================
  // Index/ID conversions

  //! Convert from indices to cell ID
  template<long D>
  CellID indices_to_cellid(const IntPoint<D>& indices) const;

  //! Convert from ID to indices
  template<long D>
  IntPoint<D> cellid_to_indices(const CellID id) const;


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
Direct_Product_Mesh_Wrapper::Direct_Product_Mesh_Wrapper(
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
// Get dimensionality of the mesh.
int Direct_Product_Mesh_Wrapper::space_dimension () const {
  return mesh_.space_dimension();
}

// ____________________________________________________________________________
// Get global mesh bounds.
template<long D>
void Direct_Product_Mesh_Wrapper::get_global_bounds(
    Point<D> *plo, Point<D> *phi) const {
  assert(D == mesh_.space_dimension());
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
// Get iterator for axis edge coordinates (beginning of array).
counting_iterator Direct_Product_Mesh_Wrapper::axis_point_begin(
    const int dim) const {
  int start_index = 0;
  return make_counting_iterator(start_index);
}

// ____________________________________________________________________________
// Get iterator for axis edge coordinates (end of array).
counting_iterator Direct_Product_Mesh_Wrapper::axis_point_end(
    const int dim) const {
  int start_index = 0;
  return make_counting_iterator(start_index + mesh_.axis_num_points(dim));
}

// ____________________________________________________________________________
// Get edge coordinate value.
double Direct_Product_Mesh_Wrapper::axis_point_coordinate(
    const int dim, const int pointid) const {
  return mesh_.axis_point_coordinate(dim, pointid);
}

// ____________________________________________________________________________
// Get number of cells along axis.
int Direct_Product_Mesh_Wrapper::axis_num_cells(const int dim) const {
  return mesh_.axis_num_points(dim) - 1;
}

// ____________________________________________________________________________
// Get number of cells in entire mesh.
int Direct_Product_Mesh_Wrapper::total_num_cells(const int dim) const {
  int count = 1;
  for (int dim = 0; dim < mesh_.space_dimension(); ++dim) {
    count *= mesh_.axis_num_points(dim) - 1;
  }
  return count;
}

// ============================================================================
// Index/ID conversions

// ____________________________________________________________________________
// Convert from indices to ID
template<long D>
CellID Direct_Product_Mesh_Wrapper::indices_to_cellid(
    const IntPoint<D>& indices) const {
  assert(D == mesh_.space_dimension());
  switch(D) {
    case 1 :
      return ((CellID) i);
      break;
    case 2 :
      CellID imax = (CellID) axis_num_cells(0);
      return ((CellID) j) * imax + ((CellID) i);
      break;
    case 3 :
      CellID imax = (CellID) axis_num_cells(0);
      CellID jmax = (CellID) axis_num_cells(1);
      return (((CellID) k) * jmax + ((CellID) j)) * imax + ((CellID) i);
      break;
  }
}

// ____________________________________________________________________________
// Convert from ID to indices
template<long D>
IntPoint<D> Direct_Product_Mesh_Wrapper::cellid_to_indices(const CellID id) const {
  assert(D == mesh_.space_dimension());
  IntPoint<D> idx;
  switch(D) {
    case 1 :
      idx[0] = (int) id;
      break;
    case 2 :
      CellID imax = (CellID) axis_num_cells(0);
      CellID i = id % imax;
      CellID j = (id - i) / imax;
      idx[0] = (int) i;
      idx[1] = (int) j;
      break;
    case 3 :
      CellID imax = (CellID) axis_num_cells(0);
      CellID jmax = (CellID) axis_num_cells(1);
      CellID i = id % imax;
      CellID j = ((id - i) / imax) % jmax;
      CellID k = (((id - i) / imax) - j) / jmax;
      idx[0] = (int) i;
      idx[1] = (int) j;
      idx[2] = (int) k;
      break;
  }
  return std::move(idx);
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
