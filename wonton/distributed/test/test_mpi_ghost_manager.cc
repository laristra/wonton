/*
 * This file is part of the Ristra portage project.
 * Please see the license file at the root of this repository, or at:
 * https://github.com/laristra/portage/blob/master/LICENSE
 */

/* NOTE: There are more comprehensive tests for this functionality in Portage */

#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "wonton/distributed/mpi_ghost_manager.h"

#include "wonton/mesh/jali/jali_mesh_wrapper.h"
#include "wonton/state/jali/jali_state_wrapper.h"

// Jali headers
#include "Mesh.hh"
#include "MeshFactory.hh"

/**
 * @brief Fixture class for ghost manager tests.
 *
 */
class GhostManagerTest : public ::testing::Test {

protected:
  using CellGhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                    Wonton::Jali_State_Wrapper,
                                                    Wonton::CELL>;
  using NodeGhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                    Wonton::Jali_State_Wrapper,
                                                    Wonton::NODE>;

  /**
   * @brief Prepare the tests.
   *
   */
  GhostManagerTest() {
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &num_ranks);
  }

  /**
   * @brief Verify send/take communication matrices.
   *
   * @param ghost_manager: the ghost manager.
   * @param mesh: the mesh.
   */
  template<class Mesh, typename GhostManager>
  void verify_communication_matrices(GhostManager const& ghost_manager,
                                     Mesh const& mesh) const {

    std::vector<std::vector<int>> send[num_ranks]; /* owned cells that are ghost for some rank */
    std::vector<std::vector<int>> take[num_ranks]; /* ghost cells that are owned by some rank */
    std::vector<int> num_sent[num_ranks];          /* number of owned cells sent by each rank */
    std::vector<int> num_take[num_ranks];          /* number of ghost cells received from each rank */
    std::vector<MPI_Request> requests;             /* list of asynchronous MPI requests */

//    GhostManager ghost_manager(mesh, state, comm);
    send[rank] = ghost_manager.owned_matrix();
    take[rank] = ghost_manager.ghost_matrix();

    for (int i = 0; i < num_ranks; ++i) {
      num_sent[rank].emplace_back(send[rank][i].size());
      num_take[rank].emplace_back(take[rank][i].size());
      // store GID for owned cells to ease verification
      for (int& id : send[rank][i]) {
        id = mesh.get_global_id(id, Wonton::CELL);
      }
    }

    if (rank > 0) {
      MPI_Request request;
      MPI_Isend(num_sent[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
      MPI_Isend(num_take[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
    } else {
      for (int i = 1; i < num_ranks; ++i) {
        MPI_Request request;
        num_sent[i].resize(num_ranks, 0);
        MPI_Irecv(num_sent[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
        num_take[i].resize(num_ranks, 0);
        MPI_Irecv(num_take[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    requests.clear();

    if (rank > 0) {
      int const offset = num_ranks * num_ranks;
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          int tag = num_ranks * (rank - 1) + i;
          MPI_Request request;
          MPI_Isend(send[rank][i].data(), send[rank][i].size(), MPI_INT, 0, tag, comm, &request);
          requests.emplace_back(request);
          MPI_Isend(take[rank][i].data(), take[rank][i].size(), MPI_INT, 0, tag + offset, comm, &request);
          requests.emplace_back(request);
        }
      }
    } else {
      // check entities count
      for (int i = 0; i < num_ranks; ++i) {
        for (int j = 0; j < num_ranks; ++j) {
          ASSERT_EQ(num_sent[i][j], num_take[j][i]);
        }
      }

      int const offset = num_ranks * num_ranks;
      for (int i = 1; i < num_ranks; ++i) {
        send[i].resize(num_ranks);
        take[i].resize(num_ranks);
        for (int j = 0; j < num_ranks; ++j) {
          if (i != j) {
            int tag = num_ranks * (i - 1) + j;
            MPI_Request request;
            send[i][j].resize(num_sent[i][j]);
            MPI_Irecv(send[i][j].data(), num_sent[i][j], MPI_INT, i, tag, comm, &request);
            requests.emplace_back(request);
            take[i][j].resize(num_take[i][j]);
            MPI_Irecv(take[i][j].data(), num_take[i][j], MPI_INT, i, tag + offset, comm, &request);
            requests.emplace_back(request);
          }
        }
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    // check that we have a perfect matching on
    // sent and received cells for each rank pair.
    if (rank == 0) {
      for (int i = 0; i < num_ranks; ++i) {
        for (int j = 0; j < num_ranks; ++j) {
          if (i == j) { continue; }
          int const num_owned = send[i][j].size();
          int const num_ghost = take[j][i].size();
          ASSERT_EQ(num_owned, num_ghost);
          for (int k = 0; k < num_owned; ++k) {
            ASSERT_EQ(send[i][j][k], take[j][i][k]);
          }
        }
      }
    }

    MPI_Barrier(comm);
  }

  /**
   * @brief Verify updated ghost values.
   *
   * @param ghost_manager: the ghost manager
   * @param rank: the current MPI rank.
   * @param num_ranks: the number of MPI ranks.
   * @param comm: the MPI communicator.
   */
  template<typename GhostManager>
  void verify_values(GhostManager const& ghost_manager) const {

    // gather all sent/received values from all ranks to the master
    std::vector<std::vector<double>> send[num_ranks]; /* values of owned cells that are ghost for some rank */
    std::vector<std::vector<double>> take[num_ranks]; /* values of ghost cells that are owned by some rank */
    std::vector<int> num_sent[num_ranks];             /* number of owned cells sent by each rank */
    std::vector<int> num_take[num_ranks];             /* number of ghost cells receive from each rank */
    std::vector<MPI_Request> requests;                /* list of asynchronous MPI requests */

    send[rank] = ghost_manager.owned_values();
    take[rank] = ghost_manager.ghost_values();

    for (int i = 0; i < num_ranks; ++i) {
      num_sent[rank].emplace_back(send[rank][i].size());
      num_take[rank].emplace_back(take[rank][i].size());
    }

    if (rank > 0) {
      MPI_Request request;
      MPI_Isend(num_sent[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
      MPI_Isend(num_take[rank].data(), num_ranks, MPI_INT, 0, rank, comm, &request);
      requests.emplace_back(request);
    } else {
      for (int i = 1; i < num_ranks; ++i) {
        MPI_Request request;
        num_sent[i].resize(num_ranks, 0);
        MPI_Irecv(num_sent[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
        num_take[i].resize(num_ranks, 0);
        MPI_Irecv(num_take[i].data(), num_ranks, MPI_INT, i, i, comm, &request);
        requests.emplace_back(request);
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);
    requests.clear();

    if (rank > 0) {
      int const offset = num_ranks * num_ranks;
      for (int i = 0; i < num_ranks; ++i) {
        if (i != rank) {
          int tag = num_ranks * (rank - 1) + i;
          MPI_Request request;
          MPI_Isend(send[rank][i].data(), send[rank][i].size(), MPI_DOUBLE, 0, tag, comm, &request);
          requests.emplace_back(request);
          MPI_Isend(take[rank][i].data(), take[rank][i].size(), MPI_DOUBLE, 0, tag + offset, comm, &request);
          requests.emplace_back(request);
        }
      }
    } else {
      int const offset = num_ranks * num_ranks;
      for (int i = 1; i < num_ranks; ++i) {
        send[i].resize(num_ranks);
        take[i].resize(num_ranks);
        for (int j = 0; j < num_ranks; ++j) {
          if (i != j) {
            int tag = num_ranks * (i - 1) + j;
            MPI_Request request;
            send[i][j].resize(num_sent[i][j]);
            MPI_Irecv(send[i][j].data(), num_sent[i][j], MPI_DOUBLE, i, tag, comm, &request);
            requests.emplace_back(request);
            take[i][j].resize(num_take[i][j]);
            MPI_Irecv(take[i][j].data(), num_take[i][j], MPI_DOUBLE, i, tag + offset, comm, &request);
            requests.emplace_back(request);
          }
        }
      }
    }

    MPI_Waitall(requests.size(), requests.data(), MPI_STATUS_IGNORE);

    if (rank == 0) {
      for (int i = 0; i < num_ranks; ++i) {
        for (int j = 0; j < num_ranks; ++j) {
          if (i == j) { continue; }
          int const num_owned = send[i][j].size();
          int const num_ghost = take[j][i].size();
          ASSERT_EQ(num_owned, num_ghost);
          for (int k = 0; k < num_owned; ++k) {
            ASSERT_DOUBLE_EQ(send[i][j][k], take[j][i][k]);
          }
        }
      }
    }

    MPI_Barrier(comm);
  }

  /** current MPI rank */
  int rank = 0;
  /** number of MPI ranks */
  int num_ranks = 1;
  /** MPI communicator */
  MPI_Comm comm = MPI_COMM_WORLD;
  /** number of owned cells */
  int num_owned_cells = 0;
  /** number of ghost cells */
  int num_ghost_cells = 0;
};

TEST_F(GhostManagerTest, CommMatrices) {

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  CellGhostManager ghost_manager(mesh, state, comm);

  verify_communication_matrices(ghost_manager, mesh);
}

TEST_F(GhostManagerTest, UpdateScalarField) {

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper  mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  int const num_owned_cells = mesh.num_owned_cells();
  int const num_ghost_cells = mesh.num_ghost_cells();
  double field[num_owned_cells + num_ghost_cells];

  for (int i = 0; i < num_owned_cells; ++i) {
    Wonton::Point<2> centroid;
    mesh.cell_centroid(i, &centroid);
    field[i] = centroid[0] + 2 * centroid[1];
  }

  state.mesh_add_data<double>(Wonton::CELL, "temperature", field);

  CellGhostManager ghost_manager(mesh, state, comm);
  ghost_manager.update_values("temperature", true);

  verify_values(ghost_manager);
}

TEST_F(GhostManagerTest, UpdateVectorField) {

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper  mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  int const num_owned_nodes = mesh.num_owned_nodes();
  int const num_ghost_nodes = mesh.num_ghost_nodes();
  Wonton::Vector<2> field[num_owned_nodes + num_ghost_nodes];

  for (int i = 0; i < num_owned_nodes; ++i) {
    Wonton::Point<2> coord;
    mesh.node_get_coordinates(i, &coord);
    field[i][0] = coord[0] + 2 * coord[1];
    field[i][1] = coord[1] + 2 * coord[0];
  }

  state.mesh_add_data<Wonton::Vector<2>>(Wonton::NODE, "velocity", field);
  
  NodeGhostManager ghost_manager(mesh, state, comm);
  ghost_manager.update_values("velocity", true);

  verify_values(ghost_manager);
}
