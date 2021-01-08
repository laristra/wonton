/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef PORTAGE_FIXED_D_MESH_WRAPPER_H_
#define PORTAGE_FIXED_D_MESH_WRAPPER_H_

#include <vector>
#include <algorithm>

#include "wonton/mesh/AuxMeshTopology.h"
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

/*
* This is compile-only test. It verifies that a fixed dimensional mesh
* wrapper is supported by the code. The actual return numbers do not
* matter but corresponds to a 2x2x2 cubic mesh.
*/

namespace Wonton {

class FixedD_Mesh_Wrapper : public AuxMeshTopology<FixedD_Mesh_Wrapper> {
 public:
  explicit FixedD_Mesh_Wrapper()
    : AuxMeshTopology<FixedD_Mesh_Wrapper>(true, true, true) {
      AuxMeshTopology<FixedD_Mesh_Wrapper>::build_aux_entities();
  }
  ~FixedD_Mesh_Wrapper() = default;

  // Requid mesh API 
  // -- counts 
  int space_dimension() const { return 3; }

  int num_owned_cells() const { return 8; }
  int num_owned_faces() const { return 36; }
  int num_owned_sides() const { return 192; }
  int num_owned_nodes() const { return 27; }
  int num_owned_corners() const { return 64; }

  int num_ghost_cells() const { return 0; }
  int num_ghost_faces() const { return 0; }
  int num_ghost_sides() const { return 0; }
  int num_ghost_nodes() const { return 0; }

  Entity_type cell_get_type(int c) const { return Entity_type::PARALLEL_OWNED; }
  Entity_type node_get_type(int n) const { return Entity_type::PARALLEL_OWNED; }
  Element_type cell_get_element_type(int c) const { return Element_type::HEX; }

  // -- connectivity information
  void cell_get_faces_and_dirs(int c, std::vector<int>* faces, std::vector<int>* dirs) const {
    faces->resize(6, 0);
    dirs->resize(6, 0);
  }
  void cell_get_nodes(int c, std::vector<int>* nodes) const { nodes->resize(8, 0); }
  void face_get_nodes(int f, std::vector<int>* nodes) const { nodes->resize(4, 0); }
  void node_get_cells(int n, const Entity_type ptype, std::vector<int>* cells) const { cells->resize(8, 0); }

  GID_t get_global_id(int id, const Entity_kind kind) const { return id; }

  // -- geometry
  void cell_get_coordinates(int c, std::vector<Point<3>> *ccoords) const { ccoords->resize(8, 0.0); }

  // -- functions that are expected to work for all dimensions
  template <int D>
  void node_get_coordinates(int n, Point<D>* pp) const { 
    if (D != 3) assert(false);
    *pp = Wonton::Point<D>();
  }
}; 

}  // namespace Wonton

#endif
