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

TEST(GhostManager, CommMatrices) {

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);

  ASSERT_GE(num_ranks, 1);

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 7, 6);
  auto jali_state = Jali::State::create(jali_mesh);
  Wonton::Jali_Mesh_Wrapper mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  using GhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                 Wonton::Jali_State_Wrapper,
                                                 Wonton::CELL>;

  std::vector<std::vector<int>> send[num_ranks]; /* owned cells that are ghost for some rank */
  std::vector<std::vector<int>> take[num_ranks]; /* ghost cells that are owned by some rank */
  std::vector<int> num_sent[num_ranks];          /* number of owned cells sent by each rank */
  std::vector<int> num_take[num_ranks];          /* number of ghost cells received from each rank */
  std::vector<MPI_Request> requests;             /* list of asynchronous MPI requests */

  GhostManager ghost_manager(mesh, state, comm);
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

TEST(GhostManager, UpdateScalarField) {

  using GhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                Wonton::Jali_State_Wrapper,
                                                Wonton::CELL>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_state = Jali::State::create(jali_mesh);

  Wonton::Jali_Mesh_Wrapper  mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  int const num_owned_cells = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  int const num_ghost_cells = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_GHOST);

  double field[num_owned_cells + num_ghost_cells];

  for (int i = 0; i < num_owned_cells; ++i) {
    Wonton::Point<2> centroid;
    mesh.cell_centroid(i, &centroid);
    field[i] = centroid[0] + 2 * centroid[1];
  }

  state.mesh_add_data<double>(Wonton::CELL, "temperature", field);

  // fill ghost values on source mesh
  GhostManager ghost_manager(mesh, state, comm);
  ghost_manager.update_values("temperature", true);

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

TEST(GhostManager, UpdateVectorField) {

  using GhostManager = Wonton::MPI_GhostManager<Wonton::Jali_Mesh_Wrapper,
                                                Wonton::Jali_State_Wrapper,
                                                Wonton::CELL>;

  int rank = 0;
  int num_ranks = 1;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &num_ranks);
  Wonton::MPIExecutor_type executor(comm);

  auto jali_mesh  = Jali::MeshFactory(comm)(0.0, 0.0, 1.0, 1.0, 5, 5);
  auto jali_state = Jali::State::create(jali_mesh);

  Wonton::Jali_Mesh_Wrapper  mesh(*jali_mesh);
  Wonton::Jali_State_Wrapper state(*jali_state);

  int const num_owned_cells = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_OWNED);
  int const num_ghost_cells = mesh.num_entities(Wonton::CELL, Wonton::PARALLEL_GHOST);

  Wonton::Vector<2> field[num_owned_cells + num_ghost_cells];

  for (int i = 0; i < num_owned_cells; ++i) {
    Wonton::Point<2> centroid;
    mesh.cell_centroid(i, &centroid);
    field[i][0] = centroid[0] + 2 * centroid[1];
    field[i][1] = centroid[1] + 2 * centroid[0];
  }

  // fill ghost values on source mesh
  GhostManager ghost_manager(mesh, state, comm);
  ghost_manager.update_values(field, 0, true);

  // gather all sent/received values from all ranks to the master
  // compact storage: vector components are flattened
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