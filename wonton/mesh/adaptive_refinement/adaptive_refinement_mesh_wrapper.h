/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_ADAPTIVE_REFINEMENT_MESH_WRAPPER_H_
#define WONTON_ADAPTIVE_REFINEMENT_MESH_WRAPPER_H_

#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"

#include <vector>
#include <algorithm>

#include "wonton/support/wonton.h"
#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Point.h"

/*!
  @file adaptive_refinement_mesh_wrapper.h
  @brief Wrapper for an Adaptive_Refinement_Mesh
*/

namespace Wonton {

/*!
  @class Adaptive_Refinement_Mesh_Wrapper adaptive_refinement_mesh_wrapper.h
  @brief A thin wrapper that implements Portage-relevant methods for
         Adaptive_Refinement_Mesh.

  The methods implemented are those required by select parts of the Wonton and
  Portage code.  This will expand as the list of components that this wrapper
  supports expands.
*/

template<int D, class CoordSys=DefaultCoordSys>
class Adaptive_Refinement_Mesh_Wrapper {

 public:

  // ==========================================================================
  // Constructors

  //! Default constructor (disabled)
  Adaptive_Refinement_Mesh_Wrapper() = delete;

  /*!
    @brief Constructor for the mesh wrapper.
    @param[in] mesh The Adaptive_Refinement_Mesh we wish to wrap.
  */
  explicit Adaptive_Refinement_Mesh_Wrapper(
      Adaptive_Refinement_Mesh<D,CoordSys> const & mesh);

  //! Copy constructor (disabled).
  Adaptive_Refinement_Mesh_Wrapper(
      Adaptive_Refinement_Mesh_Wrapper<D,CoordSys> const &) = delete;

  //! Assignment operator (disabled).
  Adaptive_Refinement_Mesh_Wrapper & operator=(
      Adaptive_Refinement_Mesh_Wrapper<D,CoordSys> const &) = delete;


  // ==========================================================================
  // Accessors

  //! Get dimensionality of the mesh.
  int space_dimension() const;

  //! Get number of cells owned by this processing element.
  int num_owned_cells() const;

  //! Get number of ghost cells on this processing element.
  int num_ghost_cells() const;

  //! Get lower and upper corners of cell bounding box
  void cell_get_bounds(const int id, Point<D> *plo, Point<D> *phi) const;


 private:

  // ==========================================================================
  // Class data

  //! The mesh to wrap.
  Adaptive_Refinement_Mesh<D,CoordSys> const & mesh_;

};  // class Adaptive_Refinement_Mesh_Wrapper


// ============================================================================
// Constructors

// ____________________________________________________________________________
// constructor
template<int D, class CoordSys>
Adaptive_Refinement_Mesh_Wrapper<D,CoordSys>::Adaptive_Refinement_Mesh_Wrapper(
    Adaptive_Refinement_Mesh<D,CoordSys> const & mesh) :
    mesh_(mesh) {
}


// ============================================================================
// Accessors

// ____________________________________________________________________________
// Get dimensionality of the mesh.
template<int D, class CoordSys>
int Adaptive_Refinement_Mesh_Wrapper<D,CoordSys>::space_dimension() const {
  return mesh_.space_dimension();
}

// ____________________________________________________________________________
// Get number of cells owned by this processing element.
template<int D, class CoordSys>
int Adaptive_Refinement_Mesh_Wrapper<D,CoordSys>::num_owned_cells() const {
  return mesh_.num_cells();
}

// ____________________________________________________________________________
// Get number of ghost cells on this processing element.
template<int D, class CoordSys>
int Adaptive_Refinement_Mesh_Wrapper<D,CoordSys>::num_ghost_cells() const {
  return 0;
}

// ____________________________________________________________________________
// Get lower and upper corners of cell bounding box
template<int D, class CoordSys>
void Adaptive_Refinement_Mesh_Wrapper<D,CoordSys>::cell_get_bounds(
    const int id, Point<D> *plo, Point<D> *phi) const {
  auto bounding_box = mesh_.cell_get_bounds(id);
  for (int d = 0; d < D; ++d) {
    (*plo)[d] = bounding_box[d][LO];
    (*phi)[d] = bounding_box[d][HI];
  }
}

}  // namespace Wonton

#endif  // WONTON_ADAPTIVE_REFINEMENT_MESH_WRAPPER_H_
