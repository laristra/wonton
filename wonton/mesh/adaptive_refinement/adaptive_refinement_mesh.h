/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_ADAPTIVE_REFINEMENT_MESH_H_
#define WONTON_ADAPTIVE_REFINEMENT_MESH_H_

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/CellID.h"
#include "wonton/support/Point.h"

/*!
  @file adaptive_refinement_mesh.h
  @brief Definition of the Adaptive_Refinement_Mesh class.

  An Adaptive_Refinement_Mesh is a basic, serial, N-dimensional mesh whose
  cells are all axis-aligned boxes.  It is intended to be a simple mesh that
  mimics many of the important features of a cell-by-cell AMR (adaptive mesh
  refinement) mesh.
 */

namespace Wonton {

/*!
  @class Adaptive_Refinement_Mesh "adaptive_refinement_mesh.h"
  @brief A simple mesh that mimics a cell-by-cell AMR mesh.

  An Adaptive_Refinement_Mesh is a basic, serial, N-dimensional mesh whose
  cells are all axis-aligned boxes.  It is intended to be a simple mesh that
  mimics many of the important features of a cell-by-cell AMR (adaptive mesh
  refinement) mesh.

  Storage of the mesh is based on a list of cells and associated geometric
  quantities.

  The Adaptive_Refinement_Mesh is designed to implement only the necessary
  functionality to test certain components in Wonton and Portage.  As the scope
  of the tests expands, the scope of functionality of the
  Adaptive_Refinement_Mesh may also expand.
 */
template<long D>
class Adaptive_Refinement_Mesh {

 private:

  //! A convenient shorthand to describe a cell
  using cell_data_t = std::array<std::array<double,2>,D>;

  //! A convenient shorthand for the type of the cell data structure.
  using mesh_data_t = std::vector<cell_data_t>;

 public:

  // ==========================================================================
  // Constructors and destructors

  //! Default constructor (disabled)
  Adaptive_Refinement_Mesh() = delete;

  /*!
    @brief Constructor to create a Adaptive_Refinement_Mesh.
    @param[in] edges The cell edge coordinates along each axis.

    Specify the edge coordinates that delineate the cells of the mesh along
    each axis.
  */
  Adaptive_Refinement_Mesh(
      double (*func)(const Point<D>), const Point<D> plo, const Point<D> phi);

  //! Assignment operator (disabled).
  Adaptive_Refinement_Mesh & operator=(const Adaptive_Refinement_Mesh<D> &) = delete;

  //! Destructor
  ~Adaptive_Refinement_Mesh();


  // ==========================================================================
  // Accessors

  //! Get the dimensionality of the mesh.
  int space_dimension() const;

  //! Get the total number of cells in the mesh.
  int num_cells() const;

  /*!
    @brief Get the lower and upper bounds of the specified cell.
    @param[in] id The ID of the cell to query.
  */
  cell_data_t cell_get_bounds(const CellID id) const;


 private:

  // ==========================================================================
  // Private support methods

  //! Check to see if the current cell needs to be refined further
  bool check_for_refinement(const cell_data_t cell, const int level);

  //! Recursive procedure to perform the actual splitting of a cell by axis.
  mesh_data_t split_cell(const int d, const cell_data_t cell);

  // Refine a single cell.
  std::pair<mesh_data_t, std::vector<int>> refine_cell(
      const mesh_data_t cells, const std::vector<int> level, const int n);

  //! Build the mesh based on the mesh corners and the refinement function.
  void build_mesh();

  // ==========================================================================
  // Class data

  //! Shorthands to clarify storage and accessors
  static const int LO = 0;
  static const int HI = 1;

  /*!
   * @brief Function pointer to a method that describes the desired refinement
   *        level across the mesh.
   * @param[in] coords The coordinates of the point at which to evaluate the
   *                   mesh refinement function.
  */
  double (*refinement_level_)(const Point<D> coords);

  //! Dimensionality of the mesh
  int dimensionality_ = -1;

  //! Mesh corner coordinates
  cell_data_t mesh_corners_;

  //! Cell corner coordinates
  mesh_data_t cells_;

};  // class Adaptive_Refinement_Mesh


// ============================================================================
// Constructors and destructors

// ____________________________________________________________________________
// Constructor
template<long D>
Adaptive_Refinement_Mesh<D>::Adaptive_Refinement_Mesh(
    double (*func)(const Point<D>), const Point<D> plo, const Point<D> phi) {
  // Verify and save dimensionality
  assert(D > 0);
  assert(D <= 3); // Maybe the logic will extend beyond this, but I'm unsure.
  dimensionality_ = D;
  // Save refinement level function pointer
  refinement_level_ = func;
  // Save mesh corners
  for (int d = 0; d < D; ++d) {
    mesh_corners_[d][LO] = plo[d];
    mesh_corners_[d][HI] = phi[d];
  }
  // Build mesh
  build_mesh();
}

// ____________________________________________________________________________
// Destructor
template<long D>
Adaptive_Refinement_Mesh<D>::~Adaptive_Refinement_Mesh() {
  cells_.clear();
  refinement_level_ = NULL;
  for (int d = 0; d < D; ++d) {
    mesh_corners_[d][LO] = 0.0;
    mesh_corners_[d][HI] = 0.0;
  }
  dimensionality_ = -1;
}


// ============================================================================
// Private support methods

// ____________________________________________________________________________
// Check to see if the current cell needs to be refined further
template<long D>
bool Adaptive_Refinement_Mesh<D>::check_for_refinement(
    const Adaptive_Refinement_Mesh<D>::cell_data_t cell, const int level) {
  // Evaluate refinement function
  // -- In principle one could do something like integrate to get the
  //    average value, but for now we'll just sample the center of the
  //    cell.  This is equivalent to second order for Cartesian geometries.
  Point<D> center;
  for (int d = 0; d < D; ++d) {
    center[d] = 0.5 * (cell[d][LO] + cell[d][HI]);
  }
  double value = refinement_level_(center);
  // Some AMR grids enforce further restrictions, such as the requirement that
  // adjacent cells differ by no more than a single level of refinement, or
  // that you need to refine in larger patches rather than by cell, etc.  This
  // is a simple mesh, and by not imposing these constraints it actually makes
  // this mesh more general.
  return(value > level);
}

// ____________________________________________________________________________
// Recursive procedure to perform the actual splitting of a cell by axis.
template<long D>
Adaptive_Refinement_Mesh<D>::mesh_data_t
    Adaptive_Refinement_Mesh<D>::split_cell(
    const int d, const Adaptive_Refinement_Mesh<D>::cell_data_t cell) {
  // Splitting point
  double midpoint = 0.5 * (cell[d][LO] + cell[d][HI]);
  // Create the lower cell
  cell_data_t cell_lo;
  for (int d2 = 0; d2 < D; ++d2) {
    cell_lo[d2][LO] = cell[d2][LO];
    cell_lo[d2][HI] = cell[d2][HI];
  }
  cell_lo[d][HI] = midpoint;
  // Create the upper cell
  cell_data_t cell_hi;
  for (int d2 = 0; d2 < D; ++d2) {
    cell_hi[d2][LO] = cell[d2][LO];
    cell_hi[d2][HI] = cell[d2][HI];
  }
  cell_hi[d][LO] = midpoint;
  // Compose all the data together
  mesh_data_t newcells;
  const int N = 1<<(d+1);
  newcells.resize(N);
  if (d > 0) {
    // Recursive case: split along the next axis, then combine results
    auto tempcells = split_cell(d-1, cell_lo);
    for (int n = 0; n < N/2; ++n) {
      newcells[n] = tempcells[n];
    }
    tempcells = split_cell(d-1, cell_hi);
    for (int n = 0; n < N/2; ++n) {
      newcells[n+N/2] = tempcells[n];
    }
  } else {
    // Base case: return the two new cells
    newcells[0] = cell_lo;
    newcells[1] = cell_hi;
  }
  return(std::move(newcells));
}

// ____________________________________________________________________________
// Refine a single cell.
// -- Remove the specified cell and replace it with refined cells.
// -- Refinement is done by splitting in half along each axis.
template<long D>
std::pair<typename Adaptive_Refinement_Mesh<D>::mesh_data_t, std::vector<int>>
    Adaptive_Refinement_Mesh<D>::refine_cell(
    const Adaptive_Refinement_Mesh<D>::mesh_data_t cells,
    const std::vector<int> level, const int n) {
  assert(n >= 0);
  assert(n < level.size());
  // Create new storage
  // -- New storage replaces 1 cell with 2^D cells
  const int num_new_cells = 1<<D - 1;
  const int newsize = level.size() + num_new_cells;
  Adaptive_Refinement_Mesh<D>::mesh_data_t newcells;
  std::vector<int> newlevel;
  newcells.resize(newsize);
  newlevel.resize(newsize);
  // Copy elements prior to the one to refine
  for (int i = 0; i < n; ++i) {
    newcells[i] = cells[i];
    newlevel[i] = level[i];
  }
  // Replace specified cell with new cells
  Adaptive_Refinement_Mesh<D>::mesh_data_t tempcells = split_cell(D-1,cells[n]);
  for (int i = 0; i < tempcells.size(); ++i) {
    newcells[n+i] = tempcells[i];
    newlevel[n+i] = level[n] + 1;
  }
  // Copy elements after the one to refine
  for (int i = n+1; i < level.size(); ++i) {
    newcells[i+num_new_cells] = cells[i];
    newlevel[i+num_new_cells] = level[i];
  }
  // Return
  return(std::move(std::make_pair(newcells, newlevel)));
}

// ____________________________________________________________________________
// Build the mesh based on the mesh corners and the refinement function.
template<long D>
void Adaptive_Refinement_Mesh<D>::build_mesh() {
  // Build the level-zero grid
  cells_.resize(1);
  cells_[0] = mesh_corners_;
  std::vector<int> level = {0};
  // Refinement loop
  bool do_refinement;
  /*do {
    // Should we do another refinement?
    do_refinement = false;
    // Which cells need to be refined?
    std::vector<bool> refine_flag;
    refine_flag.resize(level.size(), false);
    for (int n = 0; n < level.size(); ++n) {
      // Should we refine this cell?
      refine_flag[n] = check_for_refinement(cells_[n], level[n]);
      do_refinement = (do_refinement || refine_flag[n]);
    }
    // Apply refinement
    if (do_refinement)
      std::tie(cells_, level) = refine1(cells_, refine_flag, level);
  } while(do_refinement);*/
  for (int n = 0; n < level.size(); ++n) {
    if (check_for_refinement(cells_[n], level[n])) {
      std::tie(cells_, level) = refine_cell(cells_, level, n);
      --n;  // Decrement to repeat the current cell, as it was replaced
    }
  }
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get the dimensionality of the mesh.
template<long D>
int Adaptive_Refinement_Mesh<D>::space_dimension() const {
  return dimensionality_;
}

// ____________________________________________________________________________
// Get the total number of cells in the mesh.
template<long D>
int Adaptive_Refinement_Mesh<D>::num_cells() const {
  return cells_.size();
}

// ____________________________________________________________________________
// Get the lower and upper bounds of the specified cell.
template<long D>
Adaptive_Refinement_Mesh<D>::cell_data_t
    Adaptive_Refinement_Mesh<D>::cell_get_bounds(const CellID id) const {
  int n = (int) id;
  return(cells_[n]);
}

}  // namespace Wonton

#endif  // WONTON_ADAPTIVE_REFINEMENT_MESH_H_
