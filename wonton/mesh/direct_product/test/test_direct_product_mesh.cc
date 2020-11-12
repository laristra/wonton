/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/direct_product/direct_product_mesh.h"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

using Wonton::PARALLEL_OWNED;
using Wonton::PARALLEL_GHOST;
using Wonton::ALL;

// ============================================================================
/*!
  @file test_direct_product_mesh.cc
  @brief Tests for the direct_product mesh definition in direct_product_mesh.h

*/

// ============================================================================


TEST(Direct_Product_Mesh, Parallel1D) {
  double Inf = std::numeric_limits<double>::infinity();

  bool distributed = false;
  
  Wonton::Executor_type *executor;
#ifdef WONTON_ENABLE_MPI
  Wonton::MPIExecutor_type parallel_executor(MPI_COMM_WORLD);
  executor = &parallel_executor;

  // Parallel checks will fail if we are not running this on 4 ranks
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ASSERT_TRUE(nprocs == 1 || nprocs == 4);
  if (nprocs > 1) distributed = true;
#else
  Wonton::SerialExecutor_type serial_executor;
  executor = &serial_executor;
#endif

  const int D = 1;
  int MaxNG = distributed ? 2 : 0;  // max number of ghost layers to test
  for (int NG = 0; NG <= MaxNG; NG++) {  // num ghost layers

    // Create a 4 cell mesh with uneven spacing
    //                                        2^-4    2^-3   2^-2  2^-1 2^0
    const std::vector<double> axis_points1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
    int npglobal = 5;
    
    std::array<std::vector<double>,D> axis_points;
    for (int d = 0; d < D; ++d)
      axis_points[d] = axis_points1;
    Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor, NG);

    // Dimensionality
    ASSERT_EQ(mesh.space_dimension(), D);

    // Point counts along each axis (since points shared by processors
    // are owned on some and ghost on others, we cannot rely on separate
    // owned, ghost point counts)

    int NPOINT_ALL = distributed ? 2+2*NG : 5;
    
    // Default query - ALL points
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NPOINT_ALL, mesh.num_axis_points(d));

    // ALL cells queried explicitly
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NPOINT_ALL, mesh.num_axis_points(d, ALL));

  
    // Cell counts along each axis
    int NCELL_OWNED = distributed ? 1 : 4;
    int NCELL_GHOST = 2*NG;
    int NCELL_ALL = NCELL_OWNED + NCELL_GHOST;
  
    // Default query - ALL
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELL_ALL, mesh.num_axis_cells(d));

    // All cells queried explicitly
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELL_ALL, mesh.num_axis_cells(d, ALL));

    // Owned cell counts
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELL_OWNED, mesh.num_axis_cells(d, PARALLEL_OWNED));

    // Ghost cell counts
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELL_GHOST, mesh.num_axis_cells(d, PARALLEL_GHOST));

    // Verify we are correctly identifying exterior points
    for (int d = 0; d < D; ++d) {
      double lo, hi;
      mesh.get_global_coord_bounds(d, &lo, &hi);

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = distributed ? n-NG : n; // ghost point starts at -NG not 0
        double coord = mesh.get_axis_point(d, id);
        if (fabs(coord) == Inf || coord == lo || coord == hi)
          ASSERT_TRUE(mesh.point_on_external_boundary(d, id));
        else
          ASSERT_FALSE(mesh.point_on_external_boundary(d, id));
      }
    }
  
    // Cell coordinates (owned coordinate index starts from 0, ghost from -NG)
    for (int d = 0; d < D; ++d) {
      double plo, phi;
      mesh.get_local_coord_bounds(d, &plo, &phi);  // excludes ghost layers

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = n-NG;
        double expval1, expval2;

        int64_t gid = mesh.point_global_index(d, id);
        if (gid < 0) {
          expval1 = -Inf;
          expval2 = -Inf;
        } else if (gid >= npglobal) {
          expval1 = Inf;
          expval2 = Inf;
        } else {
          expval1 = pow(2,id)*plo;
          expval2 = phi/pow(2,npall-2*NG-1-id);
        }
        
        double coord = mesh.get_axis_point(d, id);
        ASSERT_EQ(expval1, coord);
        ASSERT_EQ(expval2, coord);
      }
    }
  }
}


TEST(Direct_Product_Mesh, Parallel2D_OneGhostLayer) {
  double Inf = std::numeric_limits<double>::infinity();

  bool distributed = false;
  
  Wonton::Executor_type *executor;
#ifdef WONTON_ENABLE_MPI
  Wonton::MPIExecutor_type parallel_executor(MPI_COMM_WORLD);
  executor = &parallel_executor;

  // Parallel checks will fail if we are not running this on 4 ranks
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ASSERT_TRUE(nprocs == 1 || nprocs == 4);
  if (nprocs > 1) distributed = true;
#else
  Wonton::SerialExecutor_type serial_executor;
  executor = &serial_executor;
#endif

  const int D = 2;

  // Create a 4x4 mesh
  const std::vector<double> axis_points1 = {0.0, 1.0, 2.0, 3.0, 4.0};
  int npglobal = 5;
  
  std::array<std::vector<double>,D> axis_points;
  for (int d = 0; d < D; ++d)
    axis_points[d] = axis_points1;

  // cell size in any direction
  double delta = 1.0;

  int MaxNG = distributed ? 2 : 0;  // max number of ghost layers to test
  for (int NG = 0; NG <= MaxNG; NG++) {  // num ghost layers

    // This will create two owned cells per rank
    Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor, NG);

    // Dimensionality
    ASSERT_EQ(mesh.space_dimension(), D);

    // Point counts along each axis (since points shared by processors
    // are owned on some and ghost on others, we cannot rely on separate
    // owned, ghost point counts)

    int NPOINTS_ALL = distributed ? 3+2*NG : 5;
  
    // Default query - ALL points
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NPOINTS_ALL, mesh.num_axis_points(d));

    // ALL points queried explicitly
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NPOINTS_ALL, mesh.num_axis_points(d, ALL));

  
    // Cell counts along each axis
    int NCELLS_OWNED = distributed ? 2 : 4;
    int NCELLS_GHOST = 2*NG;
    int NCELLS_ALL = NCELLS_OWNED + NCELLS_GHOST;
  
    // Default query - ALL cells
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELLS_ALL, mesh.num_axis_cells(d));

    // All cells queried explicitly
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELLS_ALL, mesh.num_axis_cells(d, ALL));

    // Owned cell counts
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELLS_OWNED, mesh.num_axis_cells(d, PARALLEL_OWNED));

    // Ghost cell counts
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELLS_GHOST, mesh.num_axis_cells(d, PARALLEL_GHOST));

    // Verify we are correctly identifying exterior points
    for (int d = 0; d < D; ++d) {
      double lo, hi;
      mesh.get_global_coord_bounds(d, &lo, &hi);

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = distributed ? n-NG : n;
        double coord = mesh.get_axis_point(d, id);
        if (fabs(coord) == Inf || coord == lo || coord == hi)
          ASSERT_TRUE(mesh.point_on_external_boundary(d, id));
        else
          ASSERT_FALSE(mesh.point_on_external_boundary(d, id));
      }
    }
  
    // Cell coordinates (owned coordinate index starts from 0, ghost from -NG)
    for (int d = 0; d < D; ++d) {
      double plo, phi;
      mesh.get_local_coord_bounds(d, &plo, &phi);

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = n-NG;
        double expval1, expval2;

        int64_t gid = mesh.point_global_index(d, id);
        if (gid < 0) {
          expval1 = -Inf;
          expval2 = -Inf;
        } else if (gid >= npglobal) {
          expval1 = Inf;
          expval2 = Inf;
        } else {
          expval1 = plo + delta*id;
          expval2 = phi - delta*(npall-2*NG-1-id);
        }
      
        double coord = mesh.get_axis_point(d, id);
        ASSERT_EQ(expval1, coord);
        ASSERT_EQ(expval2, coord);
      }
    }
  }
}


TEST(Direct_Product_Mesh, Parallel3D) {
  double Inf = std::numeric_limits<double>::infinity();

  bool distributed = false;
  
  Wonton::Executor_type *executor;
#ifdef WONTON_ENABLE_MPI
  Wonton::MPIExecutor_type parallel_executor(MPI_COMM_WORLD);
  executor = &parallel_executor;

  // Parallel checks will fail if we are not running this on 4 ranks
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  ASSERT_TRUE(nprocs == 1 || nprocs == 4);
  if (nprocs > 1) distributed = true;
#else
  Wonton::SerialExecutor_type serial_executor;
  executor = &serial_executor;
#endif

  const int D = 3;

  int MaxNG = distributed ? 2 : 0;  // max number of ghost layers to test
  for (int NG = 0; NG <= MaxNG; NG++) {  // num ghost layers

    // Create a 4x4x2 mesh
    const std::vector<double> axis_points01 = {0.0, 1.0, 2.0, 3.0, 4.0};
    int npglobal[3] = {5, 5, 3};
    std::array<std::vector<double>,D> axis_points;
    for (int d = 0; d < 2; ++d)
      axis_points[d] = axis_points01;
    axis_points[2] = std::vector<double>({-1.0, 0.0, 1.0});

    // cell size in any direction
    double delta = 1.0;

    // Partition along two directions only (if distributed). This will
    // create 8 owned cells per rank (2x2x2 blocks) per rank.
    int dpart = distributed ? 2 : 0;
    Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor, NG, dpart);

    // Dimensionality
    ASSERT_EQ(mesh.space_dimension(), D);

    // Point counts along each axis (since points shared by processors
    // are owned on some and ghost on others, we cannot rely on separate
    // owned, ghost point counts)

    // Default query - ALL points
    for (int d = 0; d < D; ++d) {
      int NPOINTS_ALL = distributed ? 3+2*NG : 5;
      if (d == 2 && !distributed)
        NPOINTS_ALL = 3;          // only 3 points and no ghost points

      ASSERT_EQ(NPOINTS_ALL, mesh.num_axis_points(d));
    }

    // ALL points queried explicitly
    for (int d = 0; d < D; ++d) {
      int NPOINTS_ALL = distributed ? 3+2*NG : 5;
      if (d == 2 && !distributed)
        NPOINTS_ALL = 3;          // only 3 points and no ghost points
    
      ASSERT_EQ(NPOINTS_ALL, mesh.num_axis_points(d, ALL));
    }

  
    // Cell counts along each axis
    // Default query - ALL cells
    for (int d = 0; d < D; ++d) {
      int NCELLS_ALL = distributed ? 2+2*NG : 4;
      if (d == 2 && !distributed)
        NCELLS_ALL = 2;
    
      ASSERT_EQ(NCELLS_ALL, mesh.num_axis_cells(d));
    }

    // All cells queried explicitly
    for (int d = 0; d < D; ++d) {
      int NCELLS_ALL = distributed ? 2+2*NG : 4;
      if (d == 2 && !distributed)
        NCELLS_ALL = 2;
    
      ASSERT_EQ(NCELLS_ALL, mesh.num_axis_cells(d, ALL));
    }

    // Owned cell counts
    for (int d = 0; d < D; ++d) {
      int NCELLS_OWNED = distributed ? 2 : 4;
      if (d == 2)
        NCELLS_OWNED = 2;
      ASSERT_EQ(NCELLS_OWNED, mesh.num_axis_cells(d, PARALLEL_OWNED));
    }

    // Ghost cell counts
    int NCELLS_GHOST = 2*NG;
    for (int d = 0; d < D; ++d)
      ASSERT_EQ(NCELLS_GHOST, mesh.num_axis_cells(d, PARALLEL_GHOST));

    // Verify we are correctly identifying exterior points
    for (int d = 0; d < D; ++d) {
      double lo, hi;
      mesh.get_global_coord_bounds(d, &lo, &hi);

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = distributed ? n-NG : n;
        double coord = mesh.get_axis_point(d, id);
        if (fabs(coord) == Inf || coord == lo || coord == hi)
          ASSERT_TRUE(mesh.point_on_external_boundary(d, id));
        else
          ASSERT_FALSE(mesh.point_on_external_boundary(d, id));
      }
    }

  
    // Cell coordinates (owned coordinate index starts from 0, ghost from -1)
    for (int d = 0; d < D; ++d) {
      double plo, phi;
      mesh.get_local_coord_bounds(d, &plo, &phi);

      int npall = mesh.num_axis_points(d, ALL);
      for (int n = 0; n < npall; ++n) {
        int id = n-NG;
        double expval1, expval2;

        int64_t gid = mesh.point_global_index(d, id);
        if (gid < 0) {
          expval1 = -Inf;
          expval2 = -Inf;
        } else if (gid >= npglobal[d]) {
          expval1 = Inf;
          expval2 = Inf;
        } else {
          expval1 = plo + delta*id;
          expval2 = phi - delta*(npall-2*NG-1-id);
        }

        double coord = mesh.get_axis_point(d, id);
        ASSERT_EQ(expval1, coord);
        ASSERT_EQ(expval2, coord);
      }
    }
  }
}

