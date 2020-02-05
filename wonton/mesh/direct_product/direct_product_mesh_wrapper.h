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
#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Point.h"

/*!
  @file direct_product_mesh_wrapper.h
  @brief Wrapper for a Direct_Product_mesh
*/

namespace Wonton {

/*!
  @class Direct_Product_Mesh_Wrapper direct_product_mesh_wrapper.h

  @brief A wrapper that implements the prescribed interface for direct
  product meshes in Portage using the Direct_Product_Mesh
  (direct_product_mesh.h)

  Unlike the Direct_Product_Mesh class that supports queries only per
  axis, this class supports some aggregate queries (e.g. number of
  cells in the mesh, cell bounds as D-dimensional points).

  This class is meant to be used with Search and Intersect classes
  operating on direct product meshes and not classes that assume an
  unstructured mesh (e.g. SearchKDTree or IntersectR2D etc.).
*/
template<int D, class CoordSys=DefaultCoordSys>
class Direct_Product_Mesh_Wrapper {

  // TODO: Add a static assert to guarantee that the mesh and mesh wrapper are
  //       using the same coordinate system?  (Same question for AR mesh and
  //       wrapper.)

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  Direct_Product_Mesh_Wrapper() = delete;

  /*!
    @brief Constructor for the mesh wrapper.
    @param[in] mesh The Direct_Product_Mesh we wish to wrap.
  */
  explicit Direct_Product_Mesh_Wrapper(
      Direct_Product_Mesh<D,CoordSys> const & mesh);

  //! Copy constructor (disabled).
  Direct_Product_Mesh_Wrapper(
      Direct_Product_Mesh_Wrapper<D,CoordSys> const &) = delete;

  //! Assignment operator (disabled).
  Direct_Product_Mesh_Wrapper & operator=(
      Direct_Product_Mesh_Wrapper<D,CoordSys> const &) = delete;

  //! Destructor
  ~Direct_Product_Mesh_Wrapper();

  // ==========================================================================
  // Accessors

  //! Get dimensionality of the mesh.
  int space_dimension () const;

  //! Is mesh distributed?
  bool distributed () const;
  
  /*!
    @brief Get global mesh bounds.

    Because a Direct_Product_Mesh is required to be an axis-aligned box, the
    global mesh bounds are specified as a Point at the lower corner and a Point
    at the upper corner (lower and upper along all axes).

    The points must have the same size as the mesh dimensionality.
  */
  void get_global_bounds(Point<D> *plo, Point<D> *phi) const;

  //! Get number of points along axis
  int axis_num_points(const int dim,
                      const Entity_type ptype = ALL) const;
  
  //! Get iterator for axis point (beginning of array).
  counting_iterator axis_point_begin(const int dim) const;

  //! Get iterator for axis point (end of array).
  counting_iterator axis_point_end(const int dim) const;

  //! Get number of nodes owned by this processing element.
  int num_owned_nodes() const;

  //! Get number of ghost nodes on this processing element.
  int num_ghost_nodes() const;

  //! Get coordinates of node
  Point<D> get_node_coordinates(const int pointid) const;
  
  //! Get axis point value.
  double get_axis_point(const int dim, const int pointid) const;

  //! Get number of cells along axis.
  int axis_num_cells(const int dim,
                     const Entity_type ptype = ALL) const;

  //! Get number of cells owned by this processing element.
  int num_owned_cells() const;

  //! Get number of ghost cells on this processing element.
  int num_ghost_cells() const;

  //! Get lower and upper corners of cell bounding box
  void cell_get_bounds(const int id, Point<D> *plo, Point<D> *phi) const;

  //! Get coordinates of cell points
  void cell_get_coordinates(const int id, std::vector<Point<D>> *ccoords) const;

  //! Iterators on mesh entity - begin
  counting_iterator begin(Entity_kind const entity,
      Entity_type const etype = ALL) const;

  //! Iterator on mesh entity - end
  counting_iterator end(Entity_kind const entity,
      Entity_type const etype = ALL) const;

  // ==========================================================================
  // Index/ID conversions

  //! Convert from indices of lo corner of cell to cell ID
  int indices_to_cellid(const std::array<int,D>& indices) const;

  //! Convert from ID to indices of lo corner of cell
  std::array<int,D> cellid_to_indices(const int id) const;

  //! Convert from indices of node to node ID
  int indices_to_nodeid(const std::array<int,D>& indices) const;

  //! Convert from node ID to indices of node
  std::array<int,D> nodeid_to_indices(const int id) const;


 private:

  // ==========================================================================
  // Class data

  //! The mesh to wrap.
  Direct_Product_Mesh<D,CoordSys> const & mesh_;

  // ==========================================================================
  // Private methods
  
  // __________________________________________________________________________
  /*! @brief Populate the d'th coordinate of cell nodes recursively
    @param[in] cellid      ID of cell
    @param[in] lowindices  Indices of lower bounding corner of cell as given
                           by cellid_to_indices
    @param[in] d           Dimension below which to populate coordinates
    @param[in] coords      Pre-sized list of coordinates (size = pow(2,D))
  */
  void populate_cell_node_coords(int cellid, std::array<int, D> lowindices,
                                 int d, std::vector<Point<D>> *coords) const;
};  // class Direct_Product_Mesh_Wrapper


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// constructor
template<int D, class CoordSys>
Direct_Product_Mesh_Wrapper<D,CoordSys>::Direct_Product_Mesh_Wrapper(
    Direct_Product_Mesh<D,CoordSys> const & mesh) :
    mesh_(mesh) {
}

// ____________________________________________________________________________
// destructor
template<int D, class CoordSys>
Direct_Product_Mesh_Wrapper<D,CoordSys>::~Direct_Product_Mesh_Wrapper() {
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get dimensionality of the mesh.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::space_dimension() const {
  return mesh_.space_dimension();
}

// ____________________________________________________________________________
// Is the mesh distributed?
template<int D, class CoordSys>
bool Direct_Product_Mesh_Wrapper<D,CoordSys>::distributed() const {
  return mesh_.distributed();
}

// ____________________________________________________________________________
// Get global mesh bounds.
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::get_global_bounds(
    Point<D> *plo, Point<D> *phi) const {
  assert(D == mesh_.space_dimension());
  mesh_.get_global_bounds(plo, phi);
}

// ____________________________________________________________________________
// Get number of points along axis.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::axis_num_points(
    const int dim, const Entity_type ptype) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.axis_num_points(dim, ptype);
}

// ____________________________________________________________________________
// Get iterator for axis points (beginning of array).
template<int D, class CoordSys>
counting_iterator Direct_Product_Mesh_Wrapper<D,CoordSys>::axis_point_begin(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  int start_index = mesh_.distributed() ? -1 : 0;
  return make_counting_iterator(start_index);
}

// ____________________________________________________________________________
// Get iterator for axis points (end of array).
template<int D, class CoordSys>
counting_iterator Direct_Product_Mesh_Wrapper<D,CoordSys>::axis_point_end(
    const int dim) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return make_counting_iterator(*axis_point_begin(dim) +
                                mesh_.axis_num_points(dim, ALL));
}

// ____________________________________________________________________________
// Get axis point value.
template<int D, class CoordSys>
double Direct_Product_Mesh_Wrapper<D,CoordSys>::get_axis_point(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.get_axis_point(dim, pointid);
}

// ____________________________________________________________________________
// Get point coordinate
template<int D, class CoordSys>
Point<D> Direct_Product_Mesh_Wrapper<D,CoordSys>::get_node_coordinates(
    const int pointid) const {
  Point<D> pnt;
  std::array<int,D> indices;
  pointid_to_indices(pointid, &indices);
  for (int d = 0; d < D; d++) pnt[d] = mesh_.get_axis_point(d, indices[d]);
}

// ____________________________________________________________________________
// Get number of cells along axis.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::axis_num_cells(
    const int dim, const Entity_type ptype) const {
  assert(dim >= 0);
  assert(dim < mesh_.space_dimension());
  return mesh_.axis_num_cells(dim, ptype);
}

// ____________________________________________________________________________
// Get number of cells owned by this processing element.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::num_owned_cells() const {
  int count = 1;
  for (int dim = 0; dim < mesh_.space_dimension(); ++dim) {
    count *= mesh_.axis_num_cells(dim, PARALLEL_OWNED);
  }
  return count;
}

// ____________________________________________________________________________
// Get number of ghost cells on this processing element.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::num_ghost_cells() const {
  if (mesh_.distributed()) {
    int num_all_cells = 1;
    for (int dim = 0; dim < mesh_.space_dimension(); ++dim)
      num_all_cells *= mesh_.axis_num_cells(dim, ALL);
    return num_all_cells - num_owned_cells();
  } else
    return 0;
}

// ____________________________________________________________________________
// Get lower and upper corners of cell bounding box
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_get_bounds(
    const int id, Point<D> *plo, Point<D> *phi) const {
  std::array<int,D> indices = cellid_to_indices(id);
  // cell N is bounded by axis points N and N+1.
  for (int d = 0; d < D; ++d) {
    (*plo)[d] = mesh_.get_axis_point(d, indices[d]);
    (*phi)[d] = mesh_.get_axis_point(d, indices[d]+1);
  }
}

// ____________________________________________________________________________
// Get number of nodes owned by this processing element.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::num_owned_nodes() const {
  int count = 1;
  for (int dim = 0; dim < mesh_.space_dimension(); ++dim)
    count *= mesh_.axis_num_points(dim, PARALLEL_OWNED);
  return count;
}

// ____________________________________________________________________________
// Get number of ghost nodes on this processing element.
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::num_ghost_nodes() const {
  if (mesh_.distributed()) {
    int num_all_nodes = 1;
    for (int dim = 0; dim < mesh_.space_dimension(); ++dim)
      num_all_nodes *= mesh_.axis_num_points(dim, ALL);
    return num_all_nodes - num_owned_nodes();
  } else
    return 0;
}


// ____________________________________________________________________________
// Populate the d'th coordinate of cell nodes and recurse
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D, CoordSys>::populate_cell_node_coords(
    int cellid, std::array<int, D> lowindices, int d,
    std::vector<Point<D>> *coords) const {
  int ncoords = coords->size()/2;
  double plo = mesh_.get_axis_point(d, lowindices[d]);
  double phi = mesh_.get_axis_point(d, lowindices[d]+1);
  int fac = pow(2, d);
  for (int i = 0; i < ncoords; i++) {
    coords[i][d] = plo;
    coords[i+fac][d] = phi;
  }

  if (d == 0)
    return;
  else
    populate_dth_coord_of_cell_nodes(cellid, lowindices, d-1, coords);
}

// ____________________________________________________________________________
// Get coordinates of cell nodes
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_get_coordinates(
    const int cellid, std::vector<Point<D>> *ccoords) const {

  int ncoords = 1;
  for (int d = 0; d < D; d++) ncoords *= 2;

  ccoords->resize(ncoords);
  std::array<int,D> indices;
  cellid_to_indices(cellid, &indices);
  populate_cell_node_coords(cellid, indices, D, ccoords);
}

// ____________________________________________________________________________
// Iterators on mesh entity - begin
template<int D, class CoordSys>
counting_iterator Direct_Product_Mesh_Wrapper<D,CoordSys>::begin(
    Entity_kind const entity, Entity_type const etype) const {
  assert(entity == CELL || entity == NODE);
  if (mesh_.distributed())
    assert(etype == ALL);  // Only ALL cells, nodes have a continuous numbering
  return(make_counting_iterator(0));
}

// ____________________________________________________________________________
// Iterator on mesh entity - end
template<int D, class CoordSys>
counting_iterator Direct_Product_Mesh_Wrapper<D,CoordSys>::end(
    Entity_kind const entity, Entity_type const etype) const {
  // Currently only valid for cells
  assert(entity == CELL || entity == NODE);
  if (mesh_.distributed())
    assert(etype == ALL);  // Only ALL cells, nodes have a continuous numbering
  // Return iterator
  int start_index = 0;
  if (entity == CELL)
    return(make_counting_iterator(start_index +
                                  num_owned_cells() + num_ghost_cells()));
  else if (entity == NODE)
    return(make_counting_iterator(start_index +
                                  num_owned_nodes() + num_ghost_nodes()));    
}

// ============================================================================
// Index/ID conversions

// ____________________________________________________________________________
// Convert from indices to Cell ID (sadly, under this scheme, the ID
// of an owned cell will depend on whether or not the mesh has a ghost
// layer) - Cell ID will always start from 0 even if the indices are
// say (-1,-1,-1)
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::indices_to_cellid(
    const std::array<int,D>& indices) const {
  assert(D == mesh_.space_dimension());
  int id = 0;
  // Loop over all but the last dimension
  for (int d = D-1; d > 0; --d) {
    int idx = distributed() ? indices[d]+1 : indices[d];
    int mult = axis_num_cells(d-1);
    id += idx;
    id *= mult;
  }
  // Handle the last dimension
  id += (distributed() ? indices[0]+1 : indices[0]);
  // Return
  return id;
}

// ____________________________________________________________________________
// Convert from Cell ID to indices
template<int D, class CoordSys>
std::array<int,D> Direct_Product_Mesh_Wrapper<D,CoordSys>::cellid_to_indices(
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
  if (distributed())
    for (int d = 0; d < D; ++d) indices[d]--;  // ghost indices must start at -1
  // Return
  return indices;
}

// ____________________________________________________________________________
// Convert from indices to Node ID (sadly, under this scheme, the ID
// of an owned node will depend on whether or not the mesh has a ghost
// layer). The first node ID will always be 0, even if it is a ghost
// node and its indices are (-1,-1,-1)
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::indices_to_nodeid(
    const std::array<int,D>& indices) const {
  assert(D == mesh_.space_dimension());
  int id = 0;
  // Loop over all but the last dimension
  for (int d = D-1; d > 0; --d) {
    int idx = distributed() ? indices[d]+1 : indices[d];
    int mult = mesh_.axis_num_points(d-1);
    id += idx;
    id *= mult;
  }
  // Handle the last dimension
  id += (distributed() ? indices[0]+1 : indices[0]);
  // Return
  return id;
}

// ____________________________________________________________________________
// Convert from Node ID to indices
template<int D, class CoordSys>
std::array<int,D> Direct_Product_Mesh_Wrapper<D,CoordSys>::nodeid_to_indices(
    const int id) const {
  assert(D == mesh_.space_dimension());
  std::array<int,D> indices;
  int residual = id;
  // Construct the denominators
  std::array<int,D> denom;
  denom[0] = 1;
  for (int d = 1; d < D; ++d) {
    denom[d] = denom[d-1] * mesh_.axis_num_points(d-1);
  }
  // Loop over all but the last dimension
  for (int d = D-1; d > 0; --d) {
    int index = residual / denom[d];
    residual -= index * denom[d];
    indices[d] = (int) index;
  }
  // Handle the last dimension
  indices[0] = (int) residual;
  if (distributed())
    for (int d = 0; d < D; ++d) indices[d]--;  // ghost indices must start at -1
  // Return
  return indices;
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
