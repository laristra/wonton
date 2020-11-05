/*
 * This file is part of the Ristra Wonton project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/wonton/blob/master/LICENSE
 */

#ifndef WONTON_MPI_GHOST_MANAGER_H
#define WONTON_MPI_GHOST_MANAGER_H

#include <map>
#include <numeric>
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"


#ifdef WONTON_ENABLE_MPI
namespace Wonton {

/**
 * @brief A communication manager to fill values on ghost cells.
 *
 * @tparam Mesh: the mesh type.
 * @tparam State: its state type.
 * @tparam entity: the entity kind.
 */
template<typename Mesh, typename State, Wonton::Entity_kind entity>
class MPI_GhostManager {

public:
  /**
   * @brief Create and initialize the manager.
   *
   * @param in_mesh: the given mesh.
   * @param in_state: its state.
   * @param in_comm: the MPI communicator.
   */
  MPI_GhostManager(Mesh const& in_mesh, State& in_state, MPI_Comm in_comm)
    : mesh(in_mesh),
      state(in_state),
      comm(in_comm)
  { cache_comm_matrices(); }

  /**
   * @brief Delete the manager.
   *
   */
  ~MPI_GhostManager() = default;

  /**
   * @brief Send owned entities values and receive ghost entities values.
   *
   * @param field: the field to be remapped.
   * @param cache: whether to cache values or not for this field.
   */
  void update_values(std::string const& field, bool cache = false) {

    if (cache) {
      if (send.values.empty() or take.values.empty()) {
        send.values.resize(num_mats);
        take.values.resize(num_mats);
        for (int m = 0; m < num_mats; ++m) {
          send.values[m].resize(num_ranks);
          take.values[m].resize(num_ranks);
        }
      }
    }

    auto field_type = state.field_type(entity, field);
    bool multimat = (field_type == Field_type::MULTIMATERIAL_FIELD);

    if (multimat) {
      assert(entity == Wonton::CELL);

      for (int m = 1; m < num_mats; ++m) {
        if (state.get_data_type(field) == typeid(double)) {

          double* matdata = nullptr;
          state.mat_get_celldata(field, m-1, &matdata);
          update_ghost_values_mat(matdata, m-1, cache);

        } else if (state.get_data_type(field) == typeid(Vector<2>)) {

          Vector<2>* matdata = nullptr;
          state.mat_get_celldata(field, m-1, &matdata);
          update_ghost_values_mat(matdata, m-1, cache);

        } else if (state.get_data_type(field) == typeid(Vector<3>)) {

          Vector<3>* matdata = nullptr;
          state.mat_get_celldata(field, m-1, &matdata);
          update_ghost_values_mat(matdata, m-1, cache);

        }
      }
      
    } else {  // mesh field

      if (state.get_data_type(field) == typeid(double)) {

        double *meshdata = nullptr;
        state.mesh_get_data(entity, field, &meshdata);
        update_ghost_values_mesh(meshdata, cache);

      } else if (state.get_data_type(field) == typeid(Vector<2>)) {

        Vector<2>* meshdata = nullptr;
        state.mesh_get_data(entity, field, &meshdata);
        update_ghost_values_mesh(meshdata, cache);

      } else if (state.get_data_type(field) == typeid(Vector<3>)) {

        Vector<3>* meshdata = nullptr;
        state.mesh_get_data(entity, field, &meshdata);
        update_ghost_values_mesh(meshdata, cache);

      }

    }
  }


  /** @brief Update ghosts of compact array of material cell values
   *
   * @param matdata  [IN/OUT}: compact material cell data
   * @param m        [IN]    : material ID
   * @param cache    [IN]    : whether to cache values or not for this field.
   *
   * matdata is compact array of material cell values, sized to
   * nmatcells_owned + nmatcells_ghost but with only the owned values
   * filled in. The ghost values are filled in by this routine.
   *
   * material ID is 0 offset unlike in update_values call
   */
  template<typename T>
  void update_ghost_values_mat(T *matdata, int m, bool cache=false) {
    int nallent_mesh = mesh.num_entities(Wonton::CELL, Wonton::ALL);
    std::unique_ptr<T> meshdata = std::make_unique<T>(nallent_mesh);

    mat_to_mesh_values(matdata, m, meshdata.get());
    update_values(meshdata.get(), m+1, cache);
    mesh_to_mat_values(meshdata.get(), m, matdata);
  }
  
  
  /** @brief Update ghosts of compact array of material cell values
   *
   * @param meshdata  [IN/OUT}: data on mesh entities
   * @param cache: whether to cache values or not for this field.
   *
   * Sized to nowned+nghost entities but with only owned entities filled in
   */
  template<typename T>
  void update_ghost_values_mesh(T *meshdata, bool cache=false) {
    update_values(meshdata, 0, cache);
  }
  
  
  /**
   * @brief Retrieve the list of owned entities to send to each rank.
   *
   * @param m: material index.
   * @return matrix of communication.
   */
  std::vector<std::vector<int>> const& owned_matrix(int m = 0) const { return send.matrix[m]; }

  /**
   * @brief Retrieve the list of ghost entities to receive from each rank.
   *
   * @param m: material index.
   * @return matrix of communication.
   */
  std::vector<std::vector<int>> const& ghost_matrix(int m = 0) const { return take.matrix[m]; }

  /**
   * @brief Retrieve values of owned entities sent to each rank.
   *
   * @return matrix of values.
   */
  std::vector<std::vector<double>> const& owned_values(int m = 0) const { return send.values[m]; }

  /**
   * @brief Retrieve values of ghost entities received from each rank.
   *
   * @return matrix of values.
   */
  std::vector<std::vector<double>> const& ghost_values(int m = 0) const { return take.values[m]; }

private:
  /**
   * @brief Build and store communication matrices to identify
   *        which owned entities data should be sent to each rank,
   *        and which ghost entities data should be received from
   *        each rank for each material.
   */
  void cache_comm_matrices() {

    // step 0: preprocessing
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_ranks);

    num_mats = state.num_materials() + 1;
    send.matrix.resize(num_mats);
    send.count.resize(num_mats);
    take.matrix.resize(num_mats);
    take.count.resize(num_mats);

    int const num_entities = mesh.num_entities(entity, Wonton::ALL);

    for (int i = 0; i < num_entities; ++i) {
      int const& gid = mesh.get_global_id(i, entity);
      auto& exchange = is_ghost(i) ? take : send;
      exchange.lookup[gid] = i;
    }

    for (int m = 0; m < num_mats; ++m) {

      send.matrix[m].resize(num_ranks);
      take.matrix[m].resize(num_ranks);
      send.count[m].resize(num_ranks);
      take.count[m].resize(num_ranks);

      std::vector<int> entities;
      int num_ghosts[num_ranks];
      int offsets[num_ranks];
      std::vector<int> received;
      std::vector<int> gids[num_ranks];

      // step 1: cache local ghost entities list
      if (m > 0) /* multi-material */{
        assert(entity == Wonton::CELL && "only for cell-centered fields");
        std::vector<int> cells;
        state.mat_get_cells(m - 1, &cells);
        for (int const& i : cells) {
          if (is_ghost(i)) {
            int const& gid = mesh.get_global_id(i, entity);
            entities.emplace_back(gid);
          }
        }
      } else {
        for (auto&& id : take.lookup) { entities.emplace_back(id.first); }
      }


      // step 2: gather number of ghosts for all ranks and deduce offsets
      num_ghosts[rank] = entities.size();
      MPI_Allgather(num_ghosts + rank, 1, MPI_INT, num_ghosts, 1, MPI_INT, comm);

      int const total_ghosts = std::accumulate(num_ghosts, num_ghosts + num_ranks, 0);
      received.resize(total_ghosts);

      int index = 0;
      for (int i = 0; i < num_ranks; ++i) {
        offsets[i] = index;
        index += num_ghosts[i];
      }

      // step 3: gather all ghost entities on all ranks
      MPI_Allgatherv(entities.data(), num_ghosts[rank], MPI_INT, received.data(), num_ghosts, offsets, MPI_INT, comm);

      // step 4: check received ghost cells, build send matrix and send entities count.
      int num_sent[num_ranks];
      std::vector<MPI_Request> requests;

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          int const& start = offsets[i];
          int const& extent = (i < num_ranks - 1 ? offsets[i+1] : total_ghosts);
          for (int j = start; j < extent; ++j) {
            int const& gid = received[j];
            if (send.lookup.count(gid)) {
              send.matrix[m][i].emplace_back(send.lookup[gid]);
              gids[i].emplace_back(gid);
            }
          }

          MPI_Request request;
          num_sent[i] = send.matrix[m][i].size();
          MPI_Isend(num_sent + i, 1, MPI_INT, i, 0, comm, &request);
          requests.emplace_back(request);
        }
      }

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          std::cout << "[" << rank << "] -> [" << i << "]: (";
          int const num_to_send = send.matrix[m][i].size();
          for (int j = 0; j < num_to_send; ++j) {
            std::cout << mesh.get_global_id(send.matrix[m][i][j], entity);
            if (j < num_to_send - 1) {
              std::cout << ", ";
            }
          }
          std::cout << ")" << std::endl;
        }
      }
#endif

      // step 5: receive ghost count per rank
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          MPI_Request request;
          MPI_Irecv(take.count[m].data() + i, 1, MPI_INT, i, 0, comm, &request);
          requests.emplace_back(request);
        }
      }

      // step 6: send give entities, receive expected ghosts.
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and num_sent[i]) {
          MPI_Request request;
          MPI_Isend(gids[i].data(), num_sent[i], MPI_INT, i, 1, comm, &request);
          requests.emplace_back(request);
        }
      }

      MPI_Waitall(requests.size(), requests.data(), status);
      requests.clear();

      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank and take.count[m][i]) {
          MPI_Request request;
          take.matrix[m][i].resize(take.count[m][i]);
          MPI_Irecv(take.matrix[m][i].data(), take.count[m][i], MPI_INT, i, 1, comm, &request);
          requests.emplace_back(request);
        }
      }

      MPI_Waitall(requests.size(), requests.data(), status);

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          std::cout << "[" << rank << "] <- [" << i << "]: (";
          for (int j = 0; j < take.count[m][i]; ++j) {
            std::cout << take.matrix[m][i][j];
            if (j < take.count[m][i] - 1) {
              std::cout << ", ";
            }
          }
          std::cout << ")" << std::endl;
        }
      }
#endif

      // verification
      int total_received = 0;
      for (int i = 0; i < num_ranks; ++i) {
        total_received += take.matrix[m][i].size();
      }
      assert(total_received == num_ghosts[rank]);
      MPI_Barrier(comm);
    }
  }

  /**
   * @brief Check if the given entity is a ghost one.
   *
   * @param i: the index of the entity.
   * @return true if a ghost entity, false otherwise.
   */
  bool is_ghost(int i) const {
    switch (entity) {
      case Wonton::CELL: return mesh.cell_get_type(i) == Wonton::PARALLEL_GHOST;
      case Wonton::NODE: return mesh.node_get_type(i) == Wonton::PARALLEL_GHOST;
      default: return false;
    }
  }

  /** @brief convert a compact array of material cell values to a full
   * array of mesh cell values for that material
   *
   * @param matvals  [IN] : compact array of material cell values
   * @param im       [IN] : mat id (0 indexed, not offset as in update_values)
   * @param meshvals [OUT]: full array of mesh cell values
   *
   * matvals size  = num owned mat cells + num ghost mat cells
   * meshvals size = num owned mesh cells + num ghost mesh cells
   *
   * meshvals must be pre-allocated and correctly sized
   *
   * NOTE: THIS WOULD BE GOOD TO HAVE IN THE STATE MANAGER DIRECTLY
   */
  template<typename T>
  void mat_to_mesh_values(T *matvals, int im, T *meshvals) {
    std::vector<int> matcells;
    state.mat_get_cells(im, &matcells);  // gives ALL cells OWNED+GHOST
    int nmatcells = matcells.size();
    for (int i = 0; i < nmatcells; i++)
      meshvals[matcells[i]] = matvals[i];
  }
  
  /** @brief convert a compact array of material cell values to a full
   * array of mesh cell values for that material
   *
   * @param meshvals [IN] : full array of mesh cell values
   * @param im       [IN] : mat id (0 indexed, not offset as in update_values)
   * @param matvals  [OUT]: compact array of material cell values
   *
   * meshvals size = num owned mesh cells + num ghost mesh cells
   * matvals size  = num owned mat cells + num ghost mat cells
   *
   * meshvals must be pre-allocated and correctly sized
   *
   * NOTE: THIS WOULD BE GOOD TO HAVE IN THE STATE MANAGER DIRECTLY
   */
  template<typename T>
  void mesh_to_mat_values(T *meshvals, int im, T *matvals) {
    std::vector<int> matcells;
    state.mat_get_cells(im, &matcells);  // gives ALL cells OWNED+GHOST
    for (auto && c : matcells) {
      int j = state.cell_index_in_material(c, im);
      matvals[j] = meshvals[c];
    }
  }

  /**
   * @brief Send owned values and receive ghost values.
   *
   * @param data: field values array (size = num_owned_in_mesh + num_ghost_in_mesh)
   * @param m: material index.
   * @param cache: whether to cache values or not.
   *
   * Note that the data must not be the compact array of material cell
   * values; rather it has to be the full mesh cell values array
   *
   * NOTE: if the inner routine does not know anything about materials
   * per se and is just populating mesh sized array, then the comm
   * matrices may not need to be material aware. The only savings we
   * are getting is the size of the buffers we are communicating but
   * not the number of communications.  So, we can consider
   * communicating the full mesh sized buffer for materials as well
   * and doing the map from mesh cell values to mat cell values
   * outside the lowest level update_values call using the state
   * manager.
   */
  void update_values(double* data, int m, bool cache) {

    // skip if no material data for this rank
    if (data == nullptr) { return; }

    std::vector<double> buffer[num_ranks];
    std::vector<MPI_Request> requests;

    // step 1: send field values
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and not send.matrix[m][i].empty()) {
        buffer[i].clear();
        for (auto&& j : send.matrix[m][i]) {  // 'j' local entity index
          assert(not is_ghost(j));
          buffer[i].emplace_back(data[j]);
        }
        MPI_Request request;
        MPI_Isend(buffer[i].data(), buffer[i].size(), MPI_DOUBLE, i, 0, comm, &request);
        requests.emplace_back(request);

        if (cache) {
          send.values[m][i].resize(buffer[i].size());
          std::copy(buffer[i].begin(), buffer[i].end(), send.values[m][i].begin());
        }
      }
    }

    // step 2: receive field data
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and take.count[m][i]) {
        buffer[i].resize(take.count[m][i]);
        MPI_Request request;
        MPI_Irecv(buffer[i].data(), take.count[m][i], MPI_DOUBLE, i, 0, comm, &request);
        requests.emplace_back(request);
      }
    }

    // step 3: update state
    MPI_Waitall(requests.size(), requests.data(), status);

    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and take.count[m][i]) {
        for (int j = 0; j < take.count[m][i]; ++j) {
          int const& gid = take.matrix[m][i][j];
          int const& lid = take.lookup[gid];
          data[lid] = buffer[i][j];
        }
        if (cache) {
          take.values[m][i].resize(take.count[m][i]);
          std::copy(buffer[i].begin(), buffer[i].end(), take.values[m][i].begin());
        }
        buffer[i].clear();
      }
    }
  }

  /**
   * @brief Update vector field values on ghost cells.
   *
   * @tparam dim: number of components of each vector.
   * @param field: field Vector array (size = num_owned_in_mesh + num_ghost_in_mesh) 
   * @param m: the material index.
   * @param cache: whether to cache values or not (for tests).
   *
   * Note that the data must not be the compact array of material cell
   * values; rather it has to be the full mesh cell values array

   */
  template<int dim>
  void update_values(Vector<dim>* field, int m, bool cache = false) {

    static_assert(0 < dim and dim < 4, "invalid dimension");

    if (cache) {
      if (send.values.empty() or take.values.empty()) {
        send.values.resize(num_mats);
        take.values.resize(num_mats);
        for (int i = 0; i < num_mats; ++i) {
          send.values[i].resize(num_ranks);
          take.values[i].resize(num_ranks);
        }
      }
    }

    // skip if no material data for this rank
    if (field == nullptr) { return; }

    // create a MPI contiguous type for serialization
    MPI_Datatype MPI_Vector;
    MPI_Type_contiguous(dim, MPI_DOUBLE, &MPI_Vector);
    MPI_Type_commit(&MPI_Vector);

    std::vector<double> buffer[num_ranks]; // stride = dim
    std::vector<MPI_Request> requests;

    // step 1: send owned values
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and not send.matrix[m][i].empty()) {
        buffer[i].clear();
        send.count[m][i] = send.matrix[m][i].size();
        buffer[i].reserve(dim * send.count[m][i]);

        for (auto&& j : send.matrix[m][i]) {  // 'j' local entity index
          assert(not is_ghost(j));
          for (int d = 0; d < dim; ++d) {
            buffer[i].emplace_back(field[j][d]);
          }
        }

        MPI_Request request;
        MPI_Isend(buffer[i].data(), send.count[m][i], MPI_Vector, i, 0, comm, &request);
        requests.emplace_back(request);

        if (cache) {
          send.values[m][i].resize(buffer[i].size());
          std::copy(buffer[i].begin(), buffer[i].end(), send.values[m][i].begin());
        }
      }
    }

    // step 2: receive ghost values
    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and take.count[m][i]) {
        buffer[i].resize(dim * take.count[m][i]);
        MPI_Request request;
        MPI_Irecv(buffer[i].data(), take.count[m][i], MPI_Vector, i, 0, comm, &request);
        requests.emplace_back(request);
      }
    }

    // step 3: update vector field
    MPI_Waitall(requests.size(), requests.data(), status);

    for (int i = 0; i < num_ranks; ++i) {
      if (i != rank and take.count[m][i]) {
        for (int j = 0; j < take.count[m][i]; ++j) {
          int const& gid = take.matrix[m][i][j];
          int const& lid = take.lookup[gid];
          for (int d = 0; d < dim; ++d) {
            field[lid][d] = buffer[i][j * dim + d];
          }
        }

        if (cache) {
          take.values[m][i].resize(dim * take.count[m][i]);
          std::copy(buffer[i].begin(), buffer[i].end(), take.values[m][i].begin());
        }
        buffer[i].clear();
      }
    }
  }

  /**
   * @struct Data for MPI communications.
   *
   */
  struct Data {
    /** lookup table for local indices */
    std::map<int, int> lookup {};
    /** communication matrices per material */
    std::vector<std::vector<std::vector<int>>> matrix {};
    /** exchanged entities count per material */
    std::vector<std::vector<int>> count {};
    /** cached values per rank for tests, stride = dim for vector fields */
    std::vector<std::vector<std::vector<double>>> values {};
  };

  /** mesh instance */
  Mesh const& mesh;
  /** mesh state */
  State& state;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_NULL;
  /** MPI status */
  MPI_Status* status = MPI_STATUS_IGNORE;
  /** current MPI rank */
  int rank = 0;
  /** number of ranks */
  int num_ranks = 1;
  /** number of materials */
  int num_mats = 0;
  /** data to send */
  Data send;
  /** data to receive */
  Data take;
};

} // namespace Wonton
#endif // ifdef WONTON_ENABLE_MPI
#endif // ifndef WONTON_MPI_GHOST_MANAGER_H
