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
#include "wonton/support/structured_partitioner.h"

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

  This is meant to be an example of a simple structured mesh
  implementation that can be used by a wrapper class (see
  direct_product_mesh_wrapper.h) that implements the interface needed
  by Portage. All getters in this class are per axis.

  The scheme for labeling cells as PARALLEL_OWNED or PARALLEL_GHOST is
  simple - each rank has a contiguous set of owned cells bounded by
  ghost cells on either side in each direction. There may be 0 to N
  ghost layers on each side (lower and upper) in each direction

  LABELING CONVENTIONS:
  The scheme for labeling points as PARALLEL_OWNED or PARALLEL_GHOST
  on a parallel partition is PER AXIS coordinate. That means when
  considering a node as a whole, the parallel status of each of its
  coordinates must be taken into consideration - if it is
  PARALLEL_GHOST in any one direction, it is a PARALLEL_GHOST. The
  scheme is as follows:

  1. The outer most points of the ghost layer are PARALLEL_GHOST

  2. Any point in the "interior" of the partition in that axis is
  PARALLEL_OWNED (i.e., if all cells connected to it are
  PARALLEL_OWNED)

  3. The high boundary of owned cells of a partition (i.e. the axis
  point with a PARALLEL_OWNED cell before it and a PARALLEL_GHOST
  after it) is PARALLEL_OWNED

  4. The low boundary of owned cells of a partition (i.e. the axis
  point with a PARALLEL_GHOST cell before it and a PARALLEL_OWNED cell
  after it) is PARALLEL_GHOST unless it forms the low boundary of the
  global domain in which case it is PARALLEL_OWNED

  The labeling scheme for nodes and cells is shown in the figure below
  (g is ghost and o is owned) for 1 layer of ghosts.

            lower       upper       lower       upper       lower       upper
           GLOBAL     partition   partition   partition   partition    GLOBAL
          boundary    boundary    boundary    boundary    boundary    boundary
              |           |           |           |           |           |
axis pts  g   o   o   o   o   g   g   g   o   o   o   g   g   g   o   o   o   g 
          +---+---+---+---+---+   +---+---+---+---+---+   +---+---+---+---+---+ 
   cells    g   o   o   o   g       g   o   o   o   g       g   o   o   o   g 
 

   Indexing conventions:
   NG = num_ghost_layers
   N  = Number of global points in the domain (exclude ghost points we created)
   Local index = index of points along an axis in the current partition

   0. Both points and cells have one index per axis
   1. The local index for the lowest PARALLEL_OWNED POINT is 0
   2. The local indices for the PARALLEL_GHOST POINTS before point 0 are -1,
      -2, ..., -NG
   3. The local index for the lowest PARALLEL_OWNED CELL is 0
   4. The local indices for the PARALLEL_GHOST CELLS before cell 0 are -1,
      -2, ..., -NG

   5. The global index for the first globally owned POINT is 0
   6. The global index for the last globallly owned POINT is N (where N is
      the total number of global points in the input)
   6. The global indices for the ghost points before the first owned point are
      -1, -2, ..., -NG for the ghost points after the last owned point are
      N, N+1, ..., N+NG-1
   7. We don't bother with global cell indices along each axis (not needed so
   far)

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
    @param[in] axis_points       Global cell boundary coords along each axis.
    @param[in] executor          Pointer to serial or parallel executor
    @param[in] num_ghost_layers  Number of ghost layers requested in parallel
    @param[in] dpart             Number of directions to partition in parallel

    Specify the axis point coordinates that delineate the cells of the mesh
    along each axis.
  */
  Direct_Product_Mesh(const std::array<std::vector<double>,D> &axis_points_in,
                      const Wonton::Executor_type *executor = nullptr,
                      const int num_ghost_layers = 0,
                      const int dpart = D);

  //! Assignment operator (disabled).
  Direct_Product_Mesh & operator=(
      const Direct_Product_Mesh<D,CoordSys> &) = delete;

  // ==========================================================================
  // Accessors

  /*!
    @brief Get the dimensionality of the mesh.
  */
  int space_dimension() const;

  /*!
    @brief Is the mesh distributed
  */
  bool distributed() const;
  
  /*!
    @brief Get the number of axis points (cell boundary coordinates)
    along a specified axis.
    @param[in] dim     The axis to query.
    @param[in] etype   The type of points to count (See Wonton::Entity_type)
    @returns           The number of points along specified axis
  */
  int axis_num_points(const int dim,
                      const Entity_type etype = Entity_type::ALL) const;

  /*!
    @brief Number of ghost layers on any end
    @returns Number of ghost cell layers

    Note that if this returns a value of 2, it means there are 2 ghost
    cells on each end of the mesh in that direction for a total of 4. 
  */
  int num_ghost_layers() const;

  /*!
    @brief Entity type of axis point
    @param[in] dim     The axis to query
    @param[in] pointid The point to query
    @returns   Wonton entity type of point (along this axis)
  */
  Entity_type axis_point_type(const int dim, const int pointid) const;
  
  /*!
    @brief Global index of point along a specified axis
    @param[in] dim     The axis to query.
    @param[in] pointid The point to query (pointid of left ghost along
                       axis is -1)
    @returns           global index of point along axis
  */
  GID_t point_global_index(const int dim, const int pointid) const;
  
  /*!
    @brief Get the number of cell along a specified axis.  
    @param[in] dim     The axis to query.
    @param[in] etype   The type of cells to count (See Wonton::Entity_type)
    @returns           The number of cells along specified axis
  */
  int axis_num_cells(const int dim,
                     const Entity_type etype = Entity_type::ALL) const;

  /*!
    @brief Get the specified axis point (coordinate value of cell corner)
    @param[in] dim     The axis to query.

    @param[in] pointid The point to query (pointid of left ghost along
    axis is -1, of first owned point is 0 and of right ghost is n
    where n is the number of owned points along the axis

    @returns coordinate value of point along axis. If the point is a
    ghost but does not exist (e.g pointid=-1 on the left boundary),
    it will throw an exception
  */
  double get_axis_point(const int dim, const int pointid) const;

  /*!
    @brief is point on external (global) boundary along given axis

    Both the owned points on the global boundary AND their boundary
    ghost counterparts will return true in this call
  */
  bool point_on_external_boundary(const int dim, const int pointid) const;

  /*!
    @brief get global coordinate bounds along axis
    @param[in]  dim  axis along which to get bounds
    @param[out] lo   lower coordinate bound along axis
    @param[out] hi   upper coordinate bound along axis
  */
  void get_global_coord_bounds(const int dim, double *lo, double *hi) const;
    
  /*!
    @brief coordinate bounds of OWNED cells in local domain (ghost layers EXCLUDED)
    @param[in]  dim  axis along which to get bounds
    @param[out] lo   lower coordinate bound along axis
    @param[out] hi   upper coordinate bound along axis
  */
  void get_local_coord_bounds(const int dim, double *lo, double *hi) const;
  
 private:

  // ==========================================================================
  // Class data

  //! coordinate bounds
  std::array<std::array<double,2>,D> coord_bounds_global_;

  //! Number of input global points in domain along each axis
  std::array<GID_t,D> num_axis_points_global_;

  //! Number of ghost layers
  int num_ghost_layers_;

  //! Global index bounds along each axis for points of __OWNED__ cells
  //! on LOCAL DOMAIN or PARTITION
  std::array<std::array<GID_t,2>,D> gid_bounds_local_;

  //! Owned+Ghost (ALL) Axis points along the D axes
  std::array<std::vector<double>,D> axis_points_;

  //! Whether the mesh is distributed
  bool distributed_=false;

#ifdef WONTON_ENABLE_MPI
  MPI_Comm mycomm_;
  int nprocs_=1;
  int rank_=0;
#endif
  
};  // class Direct_Product_Mesh


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor

template<int D, class CoordSys>
Direct_Product_Mesh<D,CoordSys>::Direct_Product_Mesh(
    const std::array<std::vector<double>,D> &axis_points_in,
    Wonton::Executor_type const *executor, int num_ghost_layers, int dpart) {

  for (int d = 0; d < D; d++) {
    coord_bounds_global_[d][0] = axis_points_in[d][0];
    coord_bounds_global_[d][1] = axis_points_in[d][axis_points_in[d].size()-1];
    num_axis_points_global_[d] = axis_points_in[d].size();
  }
  
#ifdef WONTON_ENABLE_MPI
  
  mycomm_ = MPI_COMM_NULL;
  auto mpiexec = dynamic_cast<MPIExecutor_type const *>(executor);
  if (mpiexec && mpiexec->mpicomm != MPI_COMM_NULL) {
    mycomm_ = mpiexec->mpicomm;
    MPI_Comm_size(mycomm_, &nprocs_);
    MPI_Comm_rank(mycomm_, &rank_);
    if (nprocs_ > 1) distributed_ = true;
  }
  if (distributed_) {
    // compute partition of the global domain on each processor

    int seed = 42;  // seed for randomization to ensure we get
                    // consistent results on all processors

    std::array<GID_t, D> ncells;
    for (int d = 0; d < D; d++) ncells[d] = axis_points_in[d].size()-1;
    auto partition_bounds =
        structured_partitioner<D>(nprocs_, ncells, dpart, seed);
    
    // integer bounds of points of OWNED cells in each direction on
    // local partition. The 'partition_bounds' are returned as
    // ((xlo,ylo,zlo),(xhi,yhi,zhi)) while bounds in this class are
    // stored as ((xlo,xhi),(ylo,yhi),(zlo,zhi))
    
    for (int d = 0; d < D; d++) {
      gid_bounds_local_[d][0] = partition_bounds[rank_][0][d];
      gid_bounds_local_[d][1] = partition_bounds[rank_][1][d];
    }
    
    // Build list of axis points including ghost points
    
    num_ghost_layers_ = num_ghost_layers;
    int NG = num_ghost_layers_;

    double inf  =  std::numeric_limits<double>::infinity();
    for (int d = 0; d < D; d++) {
      GID_t lo = gid_bounds_local_[d][0];
      GID_t hi = gid_bounds_local_[d][1];
      
      // owned + ghost points
      int npall = static_cast<int>(hi-lo+1+2*NG);
      axis_points_[d].resize(npall);

      if (lo == 0)
        for (int k = 0; k < NG; k++)
          axis_points_[d][k] = -inf;
      else
        for (int k = 0; k < NG; k++)
          axis_points_[d][k] = (lo-NG+k < 0) ? -inf :
              axis_points_[d][k] = axis_points_in[d][lo-NG+k];
        
      std::copy(axis_points_in[d].begin()+lo, axis_points_in[d].begin()+hi+1,
                axis_points_[d].begin()+NG);

      if (hi == num_axis_points_global_[d]-1)
        for (int k = 0; k < NG; k++)
          axis_points_[d][npall-NG+k] = inf;
      else
        for (int k = 0; k < NG; k++)
          axis_points_[d][npall-NG+k] = (hi+1+k >= num_axis_points_global_[d]) ?
              inf : axis_points_in[d][hi+1+k];
    }
  } else {
    num_ghost_layers_ = 0;
    for (int d = 0; d < D; ++d)
      axis_points_[d] = axis_points_in[d];
    for (int d = 0; d < D; ++d) {
      gid_bounds_local_[d][0] = 0;
      gid_bounds_local_[d][1] = num_axis_points_global_[d];
    }
  }
    
#else
  num_ghost_layers_ = 0;
  for (int d = 0; d < D; ++d)
    axis_points_[d] = axis_points_in[d];
    for (int d = 0; d < D; ++d) {
      gid_bounds_local_[d][0] = 0;
      gid_bounds_local_[d][1] = num_axis_points_global_[d];
    }
#endif

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
// Is the mesh distributed?
template<int D, class CoordSys>
bool Direct_Product_Mesh<D,CoordSys>::distributed() const {
  return(distributed_);
}

// ____________________________________________________________________________
// Get the number of axis points (cell boundary coordinates) of particular type
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::axis_num_points(const int dim,
                                                     Entity_type etype) const {
  assert(dim >= 0);
  assert(dim < D);
  switch(etype) {
    case PARALLEL_OWNED: {
      // Eliminate the guaranteed ghost points (i.e. lower points of
      // ghost cells before the owned cells and the upper points of
      // ghost cells after the owned cells). For instance, if there
      // are 2 ghost layers, then 4 points are guaranteed to be
      // ghosts as shown here (g is ghost, o is owned):
      //
      // axis pts  g     g    g/o    o     o     o     g     g
      //           *-----*-----*-----*-----*-----*-----*-----*
      // cells        g     g     o     o     o     o     o
      
      int n = axis_points_[dim].size()-2*num_ghost_layers_;
      
      // the upper coordinate of the partition in any direction is
      // ALWAYS owned by the current partition while the lower
      // coordinate of the partition in any direction is ghost (owned
      // by another partition) _unless_ it is on the outer low
      // boundary of global domain
      
      if (gid_bounds_local_[dim][0] != 0)
        n--;
      return n;
    }
    case PARALLEL_GHOST: {
      int n = 2*num_ghost_layers_;
      if (gid_bounds_local_[dim][0] != 0)  // see figure above
        n++;
      return n;
    }
    case ALL:
      return axis_points_[dim].size();
    default:
      return 0;
  }
}


// ____________________________________________________________________________
// Get the number of ghost layers on this mesh
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::num_ghost_layers() const {
  return num_ghost_layers_;
}

// ____________________________________________________________________________
// Get the specified axis point (cell boundary coordinate).
template<int D, class CoordSys>
double Direct_Product_Mesh<D,CoordSys>::get_axis_point(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  int npall = axis_points_[dim].size();
  assert(pointid >= -num_ghost_layers_);
  assert(pointid <= npall);

  // pointid starts at -num_ghost_layers for lowest ghost point  
  return axis_points_[dim][pointid+num_ghost_layers_];
}

// ____________________________________________________________________________
// is point on external (global) boundary along given axis
// Both the owned points on the global boundary AND their boundary
// ghost counterparts will return true in this call
template<int D, class CoordSys>
bool Direct_Product_Mesh<D,CoordSys>::point_on_external_boundary(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  GID_t global_point_index = point_global_index(dim, pointid);
  if (global_point_index <= 0)
    return true;               // coord on or outside of lower global bounds
  else if (global_point_index >= num_axis_points_global_[dim]-1)
    return true;               // coord on or outside of upper global bounds
  else
    return false;
}


// ____________________________________________________________________________
// Entity type of a point along the axis
template<int D, class CoordSys>
Entity_type Direct_Product_Mesh<D,CoordSys>::axis_point_type(
    const int dim, const int pointid) const {
  int npall = axis_num_points(dim, ALL);
  if (point_on_external_boundary(dim, pointid)) {  // global boundary
    if (pointid < 0 || pointid > npall-2*num_ghost_layers_-1)
      return BOUNDARY_GHOST;
    else if (pointid == 0)
      return PARALLEL_OWNED;  // lo bndry of owned cells in partition
    else if (pointid == npall-2*num_ghost_layers_-1)
      return PARALLEL_GHOST;  // hi bndry of owned cells in partition
  } else {
    if (pointid < 0 || pointid >= npall-2*num_ghost_layers_-1)
      return PARALLEL_GHOST;
    else
      return PARALLEL_OWNED;
  }
  return TYPE_UNKNOWN;
}
  
// ____________________________________________________________________________
// get global index of point along given axis

// Note: we expect the number of points along any one axis to be well
// within the range representable by an int
template<int D, class CoordSys>
GID_t Direct_Product_Mesh<D,CoordSys>::point_global_index(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  return gid_bounds_local_[dim][0] + pointid;
}

// ____________________________________________________________________________
// number of cells along given axis
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::axis_num_cells(
    const int dim, const Entity_type etype) const {
  assert(dim >= 0);
  assert(dim < D);
  switch (etype) {
    case ALL: 
      return axis_num_points(dim, ALL) - 1;
    case PARALLEL_OWNED:
      return axis_num_points(dim, ALL) - 1 - 2*num_ghost_layers_;
    case PARALLEL_GHOST:
      return 2*num_ghost_layers_;
    default:
      return 0;
  }
}


// ____________________________________________________________________________
// Coordinate bounds of global domain
template<int D, class CoordSys>
void Direct_Product_Mesh<D,CoordSys>::get_global_coord_bounds(const int dim,
                                                              double *plo,
                                                              double *phi) const {
    *plo = coord_bounds_global_[dim][0];
    *phi = coord_bounds_global_[dim][1];
}
    
// ____________________________________________________________________________
// Coordinate bounds of owned cells in local domain (ghost layers EXCLUDED)
template<int D, class CoordSys>
void Direct_Product_Mesh<D,CoordSys>::get_local_coord_bounds(const int dim,
                                                             double *plo,
                                                             double *phi) const {
  int npall = axis_points_[dim].size();
  *plo = axis_points_[dim][num_ghost_layers_];
  *phi = axis_points_[dim][npall-1-num_ghost_layers_];
}
    

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
