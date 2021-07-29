/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef CONDUIT_MESH_WRAPPER_H_
#define CONDUIT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <utility>

#include "conduit/conduit.hpp"                  // Conduit header
#include "conduit/conduit_blueprint.hpp"        // Conduit header

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/mesh/AuxMeshTopology.h"

namespace Wonton {

/**
 * @brief An example of mesh wrapper for Wonton.
 *
 * It includes the minimal set of methods for mesh queries.
 * It derives from Wonton::AuxMeshTopology which builds derived
 * mesh entities and connectivities (such as corners and wedges)
 * from the existing entities, and also helps to compute geometrical
 * quantities such as cell volumes and centroids.
 * It implements the "curiously recurring template pattern" to allow
 * the mesh wrapper itself to answer queries within Wonton::AuxMeshTopology.
 * see: https://en.m.wikipedia.org/wiki/Curiously_recurring_template_pattern
 *
 * @tparam Mesh: the type of the mesh to be wrapped.
 */
class Conduit_Mesh_Wrapper : public AuxMeshTopology<Conduit_Mesh_Wrapper> {
public:
  /**
   * @brief Create an instance.
   *
   * @param mesh: the mesh to be wrapped.
   */
  Conduit_Mesh_Wrapper(conduit::Node const& mesh) : mesh_(mesh) {}

  /**
   * @brief Delete the instance.
   */
  ~Conduit_Mesh_Wrapper() = default;

  /**
   * @brief Get the mesh dimension.
   *
   * @return the mesh dimension.
   */
  int space_dimension() const { return 0; }

  /**
   * @brief Get the number of owned cells.
   *
   * @return the number of cells owned by the current MPI rank.
   */
  int num_owned_cells() const { return 0; }

  /**
   * @brief Get the number of owned faces.
   *
   * @return the number of faces owned by the current MPI rank.
   */
  int num_owned_faces() const { return 0; }

  /**
   * @brief Get the number of owned nodes.
   *
   * @return the number of nodes owned by the current MPI rank.
   */
  int num_owned_nodes() const { return 0; }

  /**
   * @brief Get the number of ghost cells.
   *
   * @return the number of ghost cells stored by the current MPI rank.
   */
  int num_ghost_cells() const { return 0; }

  /**
   * @brief Get the number of ghost faces.
   *
   * @return the number of ghost faces stored by the current MPI rank.
   */
  int num_ghost_faces() const { return 0; }

  /**
   * @brief Get the number of ghost nodes.
   *
   * @return the number of ghost nodes stored by the current MPI rank.
   */
  int num_ghost_nodes() const { return 0; }

  /**
   * @brief Get the type of the given cell which can be owned or ghost
   *
   * @param cellid: the ID of the cell.
   * @return Entity_type::PARALLEL_OWNED|Entity_type::PARALLEL_GHOST
   */
  Entity_type cell_get_type(int const cellid) const { return PARALLEL_OWNED; }

  /**
   * @brief Get the type of the given node which can be owned or ghost
   *
   * @param nodeid: the ID of the node.
   * @return Entity_type::PARALLEL_OWNED|Entity_type::PARALLEL_GHOST
   */
  Entity_type node_get_type(int const nodeid) const { return PARALLEL_OWNED; }

  /**
   * @brief Get the element type of the given cell.
   *
   * @param cellid: the ID of the cell.
   * @return Element_type::TRI|QUAD|POLYGON|TET|PRISM|PYRAMID|HEX|POLYHEDRON
   */
  Element_type cell_get_element_type(int const cellid) const { return QUAD; }

  /**
   * @brief Get the list of faces and their normal directions for a given cell.
   *
   * @param[in] cellid: the ID of the cell.
   * @param[out] cfaces: list of ID of all faces that compose the cell.
   * @param[out] cfdirs: list of normal directions for those faces.
   */
  void cell_get_faces_and_dirs(int const cellid,
                               std::vector<int> *cfaces,
                               std::vector<int> *cfdirs) const { /* TODO */ }

  /**
   * @brief Get the nodes of the given cell.
   *
   * @param[in] cellid: the ID of the cell.
   * @param[out] nodes: the list of its nodes.
   */
  void cell_get_nodes(int const cellid, std::vector<int> *nodes) const { /* TODO */ }

  /**
   * @brief Get the nodes of the given face.
   *
   * @param[in] faceid: the ID of the face
   * @param[out] nodes: the list of its nodes.
   */
  void face_get_nodes(int const faceid, std::vector<int> *nodes) const { /* TODO */ }

  /**
   * @brief Get the list of all cells incident to a node.
   *
   * @param[in] nodeid: the ID of the node.
   * @param[in] ptype: the entity type (PARALLEL_OWNED|PARALLEL_GHOST)
   * @param[out] nodecells: the list of cells ID that are incident to the node.
   */
  void node_get_cells(int const nodeid,
                      Entity_type const ptype,
                      std::vector<int> *nodecells) const { /* TODO */ }

  /**
   * @brief Get the global ID of the given entity.
   *
   * @param id: the local ID of the entity.
   * @param kind: the entity kind NODE|FACE|CELL etc.
   * @return its global ID.
   */
  GID_t get_global_id(int const id, Entity_kind const kind) const { return 0; }

  /**
   * @brief Get the coordinates of the given node.
   *
   * @tparam dim: the dimension of the node.
   * @param[in] nodeid: the ID of the node.
   * @param[out] p: its coordinates.
   */
  template<int dim>
  void node_get_coordinates(int const nodeid, Point<dim> *p) const { /* TODO */ }

private:
  /** a reference to the mesh to be wrapped */
  conduit::Node const& mesh_;
};

} // namespace Wonton

#endif // CONDUIT_MESH_WRAPPER_H_
