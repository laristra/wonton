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
  simple - each rank has a contiguous set of owned cells bounded by a
  ghost cells on either side in each direction

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

  3. The high boundary of owned cells of a partition is PARALLEL_OWNED

  4. The low boundary of owned cells of a partition is PARALLEL_GHOST
  unless it forms the low boundary of the global domain in which case
  it is PARALLEL_OWNED

  The labeling scheme for nodes and cells is shown in the figure below

                Rank 0                               Rank 1

  g---G---o---O---o---O---o---G---g    g---G---o---O---o---O---o---G---G     

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
    @param[in] executor    Pointer to serial or parallel executor
    @param[in] dpart       Number of directions to partition in parallel

    Specify the axis point coordinates that delineate the cells of the mesh
    along each axis.
  */
  Direct_Product_Mesh(const std::array<std::vector<double>,D> &axis_points_in,
                      const Wonton::Executor_type *executor = nullptr,
                      const int dpart = D);

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
    @brief Is the mesh distributed
  */
  bool distributed() const;
  
  /*!
    @brief Get the number of points (cell boundary coordinates) along
    a specified axis.
    @param[in] dim     The axis to query.
    @param[in] etype   The type of points to count (See Wonton::Entity_type)
    @returns           The number of points along specified axis
  */
  int axis_num_points(const int dim,
                      const Entity_type etype = Entity_type::ALL) const;

  /*!
    @brief Global index of point along a specified axis
    @param[in] dim     The axis to query.
    @param[in] pointid The point to query (pointid of left ghost along
                       axis is -1)
    @returns           global index of point along axis
  */
  int point_global_index(const int dim, const int pointid) const;
  
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
    the returned value will be -Inf or +Inf (can check by comparing to
    std::numeric_limits<double>::infinity())
  */
  double get_axis_point(const int dim, const int pointid) const;

  /*!
    @brief is point on external boundary along given axis
  */
  bool point_on_external_boundary(const int dim, const int pointid) const;

  /*!
    @brief get global coordinate bounds along axis
    @param[in]  dim  axis along which to get bounds
    @param[out] lo   lower coordinate bound along axis
    @param[out] hi   upper coordinate bound along axis
  */
  void get_global_coord_bounds(const int dim,
                               double *lo, double *hi) const;
    
  /*!
    @brief get local coordinate bounds along axis (without the ghost layer)
    @param[in]  dim  axis along which to get bounds
    @param[out] lo   lower coordinate bound along axis
    @param[out] hi   upper coordinate bound along axis
  */
  void get_local_coord_bounds(const int dim,
                              double *lo, double *hi) const;
  
 private:

  // ==========================================================================
  // Class data

  //! Global coordinate bounds
  std::array<std::array<double,2>,D> global_coord_bounds_;

  //! Global index bounds
  std::array<std::array<int,2>,D> global_index_bounds_;

  //! Local index bounds
  std::array<std::array<int,2>,D> local_index_bounds_;

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
    Wonton::Executor_type const *executor, int dpart) {

  for (int d = 0; d < D; d++) {
    global_coord_bounds_[d][0] = axis_points_in[d][0];
    global_coord_bounds_[d][1] = axis_points_in[d][axis_points_in[d].size()-1];
    global_index_bounds_[d][0] = 0;
    global_index_bounds_[d][1] = axis_points_in[d].size()-1;
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

    std::array<int, D> ncells;
    for (int d = 0; d < D; d++) ncells[d] = axis_points_in[d].size()-1;
    auto partition_bounds =
        structured_partitioner<D>(nprocs_, ncells, dpart, seed);

    // integer bounds of owned cell points in each direction on local
    // partition partition_bounds are returned as
    // ((xlo,ylo,zlo),(xhi,yhi,zhi)) while bounds in this class are
    // stored as ((xlo,xhi),(ylo,yhi),(zlo,zhi))
    
    for (int d = 0; d < D; d++) {
      local_index_bounds_[d][0] = partition_bounds[rank_][0][d];
      local_index_bounds_[d][1] = partition_bounds[rank_][1][d];
    }
    
    // Build list of axis points including ghost points
    
    double inf  =  std::numeric_limits<double>::infinity();
    for (int d = 0; d < D; d++) {
      int lo = local_index_bounds_[d][0]-1;  // -1 on lower global boundary
      int hi = local_index_bounds_[d][1]+1;  // npoints_global_[d] on
                                             // upper global boundary
      int np = hi-lo+1;
      axis_points_[d].resize(np);
        
      axis_points_[d][0] = (lo == -1) ? -inf : axis_points_in[d][lo];
        
      std::copy(axis_points_in[d].begin()+lo+1, axis_points_in[d].begin()+hi,
                axis_points_[d].begin()+1);
      
      axis_points_[d][np-1] = ((hi == global_index_bounds_[d][1]+1) ?
                          inf : axis_points_in[d][hi]);
    }
  } else {
    for (int d = 0; d < D; ++d)
      axis_points_[d] = axis_points_in[d];
    local_index_bounds_ = global_index_bounds_;
  }
    
#else
    
  for (int d = 0; d < D; ++d)
    axis_points_[d] = axis_points_in[d];
  local_index_bounds_ = global_index_bounds_;
#endif

}

// ____________________________________________________________________________
// Destructor
template<int D, class CoordSys>
Direct_Product_Mesh<D,CoordSys>::~Direct_Product_Mesh() {
  for (int d = 0; d < D; ++d)
    axis_points_[d].clear();
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
// Get the dimensionality of the mesh.
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
  if (distributed_) {
    switch(etype) {
      case PARALLEL_OWNED: {
        // at least two ghosts corresponding corresponding to outer
        // points of ghost layer
        int n = axis_points_[dim].size()-2;
        
        // the upper coordinate of the partition in any direction is
        // always owned by the current partition while the lower
        // coordinate of the partition in any direction is ghost (owned
        // by another partition) _unless_ it is on the outer low
        // boundary of global domain
        if (local_index_bounds_[dim][0] != global_index_bounds_[dim][0])
          n--;
        return n;
      }
      case PARALLEL_GHOST: {
        int n = 2;
        if (local_index_bounds_[dim][0] != global_index_bounds_[dim][0])
          n++;
        return n;
      }
      case ALL:
        return axis_points_[dim].size();
      default:
        return 0;
    }
  } else {
    switch(etype) {
      case PARALLEL_OWNED: case ALL:
        return axis_points_[dim].size();
      default:
        return 0;
    }
  }
}

// ____________________________________________________________________________
// Get the specified axis point (cell boundary coordinate).
template<int D, class CoordSys>
double Direct_Product_Mesh<D,CoordSys>::get_axis_point(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  int npall = axis_points_[dim].size();
  if (distributed_) {
    assert(pointid >= -1);
    assert(pointid <= npall);
    return axis_points_[dim][pointid+1];  // starts at -1 for lower ghost point
  } else {
    assert(pointid >= 0);
    assert(pointid < npall);
    return axis_points_[dim][pointid];    // no ghosts
  }
}

// ____________________________________________________________________________
// is point external boundary along given axis
template<int D, class CoordSys>
bool Direct_Product_Mesh<D,CoordSys>::point_on_external_boundary(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  int global_point_index = point_global_index(dim, pointid);
  if (global_point_index <= global_index_bounds_[dim][0])
    return true;               // coord on or outside of lower global bounds
  else if (global_point_index >= global_index_bounds_[dim][1])
    return true;               // coord on or outside of upper global bounds
  else
    return false;
}


// ____________________________________________________________________________
// get global index of point along given axis

// Note: we expect the number of points along any one axis to be well
// within the range representable by an int
template<int D, class CoordSys>
int Direct_Product_Mesh<D,CoordSys>::point_global_index(
    const int dim, const int pointid) const {
  assert(dim >= 0);
  assert(dim < D);

  return local_index_bounds_[dim][0] + pointid;
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
      if (distributed_) {
        return axis_num_points(dim, ALL) - 3;  // ALL cells minus 2 ghosts
      } else {
        return axis_num_points(dim, ALL) - 1;  // same as ALL cells in serial
      }
    case PARALLEL_GHOST:
      return distributed_ ? 2 : 0;
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
    *plo = global_coord_bounds_[dim][0];
    *phi = global_coord_bounds_[dim][1];
}
    
// ____________________________________________________________________________
// Coordinate bounds of local domain (without the ghost layer
template<int D, class CoordSys>
void Direct_Product_Mesh<D,CoordSys>::get_local_coord_bounds(const int dim,
                                                             double *plo,
                                                             double *phi) const {
  int np = axis_points_[dim].size();
  *plo = axis_points_[dim][0];
  *phi = axis_points_[dim][np-1];
}
    

}  // namespace Wonton

#endif  // WONTON_DIRECT_PRODUCT_MESH_H_
