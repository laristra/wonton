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

  // Reference to underlying mesh
  Direct_Product_Mesh<D, CoordSys> const & mesh() const;
  
  //! Get dimensionality of the mesh.
  int space_dimension () const;

  //! Is mesh distributed?
  bool distributed () const;

  //! Number of ghost layers in each direction at each end
  int num_ghost_layers () const;
  
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

  //! Node type (PARALLEL_OWNED, PARALLEL_GHOST or BOUNDARY_GHOST)
  Entity_type node_get_type(const int id) const;  
  
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

  //! Cell type (PARALLEL_OWNED, PARALLEL_GHOST or BOUNDARY_GHOST)
  Entity_type cell_get_type(const int id) const;  
  
  //! Get lower and upper corners of cell bounding box
  void cell_get_bounds(const int id, Point<D> *plo, Point<D> *phi) const;

  //! Get coordinates of cell points
  void cell_get_coordinates(const int id, std::vector<Point<D>> *ccoords) const;

  //! Get cell volume
  double cell_volume(const int id) const;

  //! Get cell centroid
  void cell_centroid(const int id, Point<D> *ccen) const;

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


  // Adjacencies

  /*!
    @brief Get the list of cell IDs for all cells attached to a specific
    cell through its faces.
    @param[in] cellid The ID of the cell.
    @param[in] ptype The Entity_type (e.g. PARALLEL_OWNED)
    @param[out] adjcells The list of cell IDs for all cells attached
    to cell @c cellid through its faces, excluding @c cellid. The
    order of cells is consistent with the order of respective
    interfaces returned by cell_get_faces_and_dirs, but all the
    boundary faces are skipped.
  */
  void cell_get_node_adj_cells(const int cellid, const Entity_type ptype,
                               std::vector<int> *adjcells) const;

  /*!
    @brief Get the list of node IDs for all nodes attached to all cells
    attached to a specific node.
    @param[in] nodeid The ID of the node.
    @param[in] ptype The Entity_type (e.g. PARALLEL_OWNED)
    @param[out] adjnodes The list of node IDs for all cells attached to
    @c nodeid, excluding @c nodeid.
  */
  void node_get_cell_adj_nodes(const int nodeid, const Entity_type ptype,
                               std::vector<int> *adjnodes) const;

  

 private:

  // ==========================================================================
  // Class data

  //! The mesh to wrap.
  Direct_Product_Mesh<D,CoordSys> const & mesh_;

  // ==========================================================================
  // Private methods
  
  /*!  
    @brief Recursively fill in the low or high coordinates of nodes
    of a cell or its lower dimensional entities (faces/edges)

    @param[in] dim         Dimension at and below which to populate coordinates
    @param[in] indices     Indices of lower bounding corner of cell as given
                           by cellid_to_indices
    @param[in] lowval      Whether to use low bound (true) or high (false)
    @param[in] coords      Pre-sized list of coordinates

    Top level call should start with d=0
  */
  void fill_dth_coord(
      int const dim, std::array<int,D> lowindices, bool lowval,
      typename std::vector<Point<D>>::iterator coords_begin,
      typename std::vector<Point<D>>::iterator coords_end) const;

  /*!
    @brief Recursively build combinations of indices from index ranges

    @param[in] dim          Dimension at or below which to build combinations
    @param[in] indices      Partially filled array of higher dimensional indices
    @param[in] index_ranges Ranges of the indices in each dimension
    @param[in/out] indices_list List of index combinations
  */
  
  void build_index_combinations(int dim,
                                std::array<int,D> indices,
                                std::array<std::array<int,2>,D> index_ranges,
                                std::vector<std::array<int,D>> *indices_list) const;
  
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
// Reference to underlying mesh
template<int D, class CoordSys>
Direct_Product_Mesh<D, CoordSys> const &
Direct_Product_Mesh_Wrapper<D,CoordSys>::mesh() const {
  return mesh_;
}
  
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
// Number of ghost layers in each direction on each end
template<int D, class CoordSys>
int Direct_Product_Mesh_Wrapper<D,CoordSys>::num_ghost_layers() const {
  return mesh_.num_ghost_layers();
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
  int start_index = -mesh_.num_ghost_layers();  // could be 0
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
// Node type (PARALLEL_OWNED, PARALLEL_GHOST or BOUNDARY_GHOST)
template<int D, class CoordSys>
Entity_type
Direct_Product_Mesh_Wrapper<D,CoordSys>::node_get_type(const int id) const {
  std::array<int, D> indices = nodeid_to_indices(id);

  // Order of precedence is BOUNDARY_GHOST, PARALLEL_GHOST,
  // PARALLEL_OWNED i.e. if point is OWNED on one axis, PARALLEL_GHOST
  // on another axis and BOUNDARY_GHOST along a third, the point will
  // be labeled a BOUNDARY_GHOST
  Entity_type etype = PARALLEL_OWNED;  
  for (int d = 0; d < D; d++) {
    Entity_type dtype = mesh_.axis_point_type(d, indices[d]);
    if (dtype == BOUNDARY_GHOST) {
      etype = BOUNDARY_GHOST;
      break;                  // game over; this takes highest precedence
    } else if (dtype == PARALLEL_GHOST) {
      etype = PARALLEL_GHOST;
    }
  }
  return etype;
}
  
// ____________________________________________________________________________
// Cell type (PARALLEL_OWNED, PARALLEL_GHOST or BOUNDARY_GHOST)
template<int D, class CoordSys>
Entity_type
Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_get_type(const int id) const {
  std::array<int, D> indices = cellid_to_indices(id);

  // Order of precedence is BOUNDARY_GHOST, PARALLEL_GHOST,
  // PARALLEL_OWNED i.e. if point is OWNED on one axis, PARALLEL_GHOST
  // on another axis and BOUNDARY_GHOST along a third, the point will
  // be labeled a BOUNDARY_GHOST
  Entity_type etype = PARALLEL_OWNED;  
  for (int d = 0; d < D; d++) {
    int lo = indices[d];
    int hi = indices[d]+1;
    Entity_type lotype = mesh_.axis_point_type(d, lo);
    Entity_type hitype = mesh_.axis_point_type(d, lo);
    if (lotype == BOUNDARY_GHOST || hitype == BOUNDARY_GHOST) {
      etype = BOUNDARY_GHOST;  // game over; this takes highest precedence
      break;
    } else if (lotype == PARALLEL_GHOST || hitype == PARALLEL_GHOST) {
      etype = PARALLEL_GHOST;
    }
  }
  return etype;
}
  

// ____________________________________________________________________________
// Recursively fill in the low or high coordinates of nodes of a cell
// or it's lower dimensional entities (faces/edges)
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::fill_dth_coord(
    int const dim, std::array<int,D> indices, bool lowval,
    typename std::vector<Point<D>>::iterator coords_begin,
    typename std::vector<Point<D>>::iterator coords_end) const {
  int idx = lowval ? indices[dim] : indices[dim]+1;  // low or high value
  double val = mesh_.get_axis_point(dim, idx);
  for (auto it = coords_begin; it != coords_end; it++) {
    Point<D> &coord = *it;
    coord[dim] = val;
  }
  if (dim == 0)
    return;
  else {
    int ncoords = coords_end - coords_begin;
    // Fill in the low coordinates in the next dimension for a sub-entity
    fill_dth_coord(dim-1, indices, true,
                   coords_begin, coords_begin + ncoords/2);
    // Fill in the high coordinates in the next dimension for a sub-entity
    fill_dth_coord(dim-1, indices, false,
                   coords_begin + ncoords/2, coords_end);
  }
}
                                                                
// ____________________________________________________________________________
// Get coordinates of cell nodes
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_get_coordinates(
    const int cellid, std::vector<Point<D>> *ccoords) const {

  int ncoords = 1;
  for (int d = 0; d < D; d++) ncoords *= 2;

  ccoords->resize(ncoords);
  std::array<int,D> indices = cellid_to_indices(cellid);

  // Recursively fill in the low and high plane coordinates of cell
  // (true - low value, false - high value)
  fill_dth_coord(D-1, indices, true,
                 ccoords->begin(), ccoords->begin() + ncoords/2);
  fill_dth_coord(D-1, indices, false,
                 ccoords->begin() + ncoords/2, ccoords->end());
}


// ____________________________________________________________________________
// Get volume of cell accounting for coordinate system
template<int D, class CoordSys>
double Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_volume(
    const int cellid) const {
  Point<D> lo, hi;
  cell_get_bounds(cellid, &lo, &hi);

  double moment0 = 1.0;  // cartesian volume
  for (int d = 0; d < D; d++)
    moment0 *= (hi[d]-lo[d]);
  
  CoordSys::modify_volume(moment0, lo, hi);
  return moment0;
}

// ____________________________________________________________________________
// Get centroid of cell accounting for coordinate system
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_centroid(
    const int cellid, Point<D> *ccen) const {
  Point<D> lo, hi;
  cell_get_bounds(cellid, &lo, &hi);

  double moment0 = 1.0;  // cartesian volume
  for (int d = 0; d < D; d++)
    moment0 *= (hi[d]-lo[d]);

  Point<D> moment1 = 0.5*(lo+hi)*moment0;

  CoordSys::modify_volume(moment0, lo, hi);
  CoordSys::modify_first_moments(moment1, lo, hi);
  *ccen = moment1/moment0;
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
  else
    return begin(entity, etype);
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
    int idx = indices[d]+mesh_.num_ghost_layers();
    int mult = axis_num_cells(d-1);
    id += idx;
    id *= mult;
  }
  // Handle the last dimension
  id += indices[0]+mesh_.num_ghost_layers();
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
  // indices start at -NG where NG = number of ghost layers
  for (int d = 0; d < D; ++d)
    indices[d] -= mesh_.num_ghost_layers();
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
    int idx = indices[d]+mesh_.num_ghost_layers();
    int mult = mesh_.axis_num_points(d-1);
    id += idx;
    id *= mult;
  }
  // Handle the last dimension
  id += indices[0]+mesh_.num_ghost_layers();
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
  // Indices start at -NG where NG is number of ghost layers
  for (int d = 0; d < D; ++d)
    indices[d] -= mesh_.num_ghost_layers();
  // Return
  return indices;
}

// Recursively build combinations of indices from index ranges
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::build_index_combinations(
    int dim,
    std::array<int,D> indices,
    std::array<std::array<int,2>,D> index_ranges,
    std::vector<std::array<int,D>> *indices_list) const {
  for (int i = index_ranges[dim][0]; i <= index_ranges[dim][1]; i++) {
    indices[dim] = i;  // The other indices are filled in already
    if (dim == 0)
      indices_list->push_back(indices);
    else  // recurse down to the build add lower dimensional indices
      build_index_combinations(dim-1, indices, index_ranges, indices_list);
  }
}

// Get the list of cell IDs for all cells attached to a specific cell
// through its faces.
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::cell_get_node_adj_cells(
    const int cellid, const Entity_type ptype, std::vector<int> *adjcells) const {
  std::array<int, D> indices = cellid_to_indices(cellid);

  std::array<std::array<int, 2>, D> index_ranges;
  int num_ghost_layers = mesh_.num_ghost_layers();
  for (int d = 0; d < D; d++) {
    int num_points_all = mesh_.axis_num_points(d, Wonton::ALL);
    index_ranges[d][0] = (indices[d]-1 >= -num_ghost_layers) ?
        indices[d]-1 : indices[d];
    index_ranges[d][1] = (indices[d]+1 < num_points_all-num_ghost_layers-1) ?
        indices[d]+1 : indices[d];
  }

  int nmax = 1;
  for (int d = 0; d < D; d++) nmax *= 3;
  
  std::vector<std::array<int,D>> cell_indices;
  cell_indices.reserve(nmax);
  std::array<int, D> tmpindices;
  build_index_combinations(D-1, tmpindices, index_ranges, &cell_indices);

  for (auto const& inds : cell_indices) {
    int id = indices_to_cellid(inds);
    if (ptype == ALL || ptype == cell_get_type(id))
      adjcells->push_back(id);
  }
}

// Get the list of node IDs for all nodes attached to all cells
// attached to a specific node.
template<int D, class CoordSys>
void Direct_Product_Mesh_Wrapper<D,CoordSys>::node_get_cell_adj_nodes(
    const int nodeid, const Entity_type ptype, std::vector<int> *adjnodes) const {
  std::array<int, D> indices = nodeid_to_indices(nodeid);

  std::array<std::array<int, 2>, D> index_ranges;
  int num_ghost_layers = mesh_.num_ghost_layers();
  for (int d = 0; d < D; d++) {
    int num_points_all = mesh_.axis_num_points(d, Wonton::ALL);
    index_ranges[d][0] = (indices[d]-1 >= -num_ghost_layers) ?
        indices[d]-1 : indices[d];
    index_ranges[d][1] = (indices[d]+1 < num_points_all-num_ghost_layers) ?
        indices[d]+1 : indices[d];
  }

  int nmax = 1;
  for (int d = 0; d < D; d++) nmax *= 3;
  
  std::vector<std::array<int,D>> node_indices;
  node_indices.reserve(nmax);
  std::array<int, D> tmpindices;
  build_index_combinations(D-1, tmpindices, index_ranges, &node_indices);

  for (auto const& inds : node_indices) {
    int id = indices_to_nodeid(inds);
    if (ptype == ALL || ptype == node_get_type(id))
      adjnodes->push_back(id);
  }
}

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_WRAPPER_H_
