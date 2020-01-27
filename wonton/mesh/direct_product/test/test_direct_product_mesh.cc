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

#ifdef WONTON_ENABLE_MPI

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

  // Create a 4 cell mesh with uneven spacing
  //                                        2^-4    2^-3   2^-2  2^-1 2^0
  const std::vector<double> axis_points1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  std::array<std::vector<double>,D> axis_points;
  for (int d = 0; d < D; ++d)
    axis_points[d] = axis_points1;
  Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Point counts along each axis (since points shared by processors
  // are owned on some and ghost on others, we cannot rely on separate
  // owned, ghost point counts)

  int NPOINT_ALL = distributed ? 4 : 5;
    
  // Default query - ALL points
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NPOINT_ALL, mesh.axis_num_points(d));

  // ALL cells queried explicitly
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NPOINT_ALL, mesh.axis_num_points(d, ALL));

  
  // Cell counts along each axis
  int NCELL_OWNED = distributed ? 1 : 4;
  int NCELL_GHOST = distributed ? 2 : 0;
  int NCELL_ALL = NCELL_OWNED + NCELL_GHOST;
  
  // Default query - ALL
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELL_ALL, mesh.axis_num_cells(d));

  // All cells queried explicitly
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELL_ALL, mesh.axis_num_cells(d, ALL));

  // Owned cell counts
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELL_OWNED, mesh.axis_num_cells(d, PARALLEL_OWNED));

  // Ghost cell counts
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELL_GHOST, mesh.axis_num_cells(d, PARALLEL_GHOST));

  // Verify we are correctly identifying exterior points
  for (int d = 0; d < D; ++d) {
    double lo, hi;
    mesh.get_global_coord_bounds(d, &lo, &hi);

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n; // ghost point starts at -1 not 0
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

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n;
      double expval1, expval2;
      
      if (n == 0) {
        expval1 = plo;
        expval2 = plo;
      } else if (n == npall-1) {
        expval1 = phi;
        expval2 = phi;
      } else {
        expval1 = (plo != -Inf) ? pow(2,n)*plo : phi/pow(2,npall-n-1);
        expval2 = (phi !=  Inf) ? phi/pow(2,npall-n-1) : pow(2,n)*plo;
      }
        
      double coord = mesh.get_axis_point(d, id);
      ASSERT_EQ(expval1, coord);
      ASSERT_EQ(expval2, coord);
    }
  }
}


TEST(Direct_Product_Mesh, Parallel2D) {
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
  std::array<std::vector<double>,D> axis_points;
  for (int d = 0; d < D; ++d)
    axis_points[d] = axis_points1;

  // cell size in any direction
  double delta = 1.0;

  // This will create two owned cells per rank
  Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Point counts along each axis (since points shared by processors
  // are owned on some and ghost on others, we cannot rely on separate
  // owned, ghost point counts)

  int NPOINTS_ALL = distributed ? 5 : 5;  // in dist. case (3 + 2)
  
  // Default query - ALL points
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NPOINTS_ALL, mesh.axis_num_points(d));

  // ALL points queried explicitly
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NPOINTS_ALL, mesh.axis_num_points(d, ALL));

  
  // Cell counts along each axis
  int NCELLS_OWNED = distributed ? 2 : 4;
  int NCELLS_GHOST = distributed ? 2 : 0;
  int NCELLS_ALL = NCELLS_OWNED + NCELLS_GHOST;
  
  // Default query - ALL cells
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELLS_ALL, mesh.axis_num_cells(d));

  // All cells queried explicitly
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELLS_ALL, mesh.axis_num_cells(d, ALL));

  // Owned cell counts
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELLS_OWNED, mesh.axis_num_cells(d, PARALLEL_OWNED));

  // Ghost cell counts
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELLS_GHOST, mesh.axis_num_cells(d, PARALLEL_GHOST));

  // Verify we are correctly identifying exterior points
  for (int d = 0; d < D; ++d) {
    double lo, hi;
    mesh.get_global_coord_bounds(d, &lo, &hi);

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n;
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

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n;
      double expval1, expval2;

      if (n == 0) {
        expval1 = plo;
        expval2 = plo;
      } else if (n == npall-1) {
        expval1 = phi;
        expval2 = phi;
      } else {
        int global_id = mesh.point_global_index(d, id);
        int npaxis = axis_points[d].size();  // num input points along d
        expval1 = (plo != -Inf) ? plo + delta*n :
            axis_points[d][0] + delta*global_id;
        expval2 = (phi !=  Inf) ? phi - delta*(npall-1-n) :
            axis_points[d][npaxis-1] - delta*(npaxis-1-global_id);
      }
      
      double coord = mesh.get_axis_point(d, id);
      ASSERT_EQ(expval1, coord);
      ASSERT_EQ(expval2, coord);
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

  // Create a 4x4x2 mesh
  const std::vector<double> axis_points01 = {0.0, 1.0, 2.0, 3.0, 4.0};
  std::array<std::vector<double>,D> axis_points;
  for (int d = 0; d < 2; ++d)
    axis_points[d] = axis_points01;
  axis_points[2] = std::vector<double>({-1.0, 0.0, 1.0});

  // cell size in any direction
  double delta = 1.0;

  // Partition along two directions only (if distributed). This will
  // create 8 owned cells per rank (2x2x2 blocks) per rank.
  int dpart = distributed ? 2 : 0;
  Wonton::Direct_Product_Mesh<D> mesh(axis_points, executor, dpart);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Point counts along each axis (since points shared by processors
  // are owned on some and ghost on others, we cannot rely on separate
  // owned, ghost point counts)

  // Default query - ALL points
  for (int d = 0; d < D; ++d) {
    int NPOINTS_ALL = 5;  // 5 in serial case; 3 owned + 2 ghost in dist case
    if (d == 2 && !distributed)
      NPOINTS_ALL = 3;          // only 3 points and no ghost points

    ASSERT_EQ(NPOINTS_ALL, mesh.axis_num_points(d));
  }

  // ALL points queried explicitly
  for (int d = 0; d < D; ++d) {
    int NPOINTS_ALL = 5;  // 5 in serial case; 3 owned + 2 ghost in dist case
    if (d == 2 && !distributed)
      NPOINTS_ALL = 3;          // only 3 points and no ghost points
    
    ASSERT_EQ(NPOINTS_ALL, mesh.axis_num_points(d, ALL));
  }

  
  // Cell counts along each axis
  // Default query - ALL cells
  for (int d = 0; d < D; ++d) {
    int NCELLS_ALL = 4;  // 4 in serial case; 2 owned + 2 ghost in dist case
    if (d == 2 && !distributed)
      NCELLS_ALL = 2;
    
    ASSERT_EQ(NCELLS_ALL, mesh.axis_num_cells(d));
  }

  // All cells queried explicitly
  for (int d = 0; d < D; ++d) {
    int NCELLS_ALL = 4;  // 4 in serial case; 2 owned + 2 ghost in dist case
    if (d == 2 && !distributed)
      NCELLS_ALL = 2;
    
    ASSERT_EQ(NCELLS_ALL, mesh.axis_num_cells(d, ALL));
  }

  // Owned cell counts
  for (int d = 0; d < D; ++d) {
    int NCELLS_OWNED = distributed ? 2 : 4;
    if (d == 2)
      NCELLS_OWNED = 2;
    ASSERT_EQ(NCELLS_OWNED, mesh.axis_num_cells(d, PARALLEL_OWNED));
  }

  // Ghost cell counts
  int NCELLS_GHOST = distributed ? 2 : 0;
  for (int d = 0; d < D; ++d)
    ASSERT_EQ(NCELLS_GHOST, mesh.axis_num_cells(d, PARALLEL_GHOST));

  // Verify we are correctly identifying exterior points
  for (int d = 0; d < D; ++d) {
    double lo, hi;
    mesh.get_global_coord_bounds(d, &lo, &hi);

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n;
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

    int npall = mesh.axis_num_points(d, ALL);
    for (int n = 0; n < npall; ++n) {
      int id = distributed ? n-1 : n;
      double expval1, expval2;

      if (n == 0) {
        expval1 = plo;
        expval2 = plo;
      } else if (n == npall-1) {
        expval1 = phi;
        expval2 = phi;
      } else {
        int global_id = mesh.point_global_index(d, id);
        int npaxis = axis_points[d].size();  // num input points along d
        expval1 = (plo != -Inf) ? plo + delta*n :
            axis_points[d][0] + delta*global_id;
        expval2 = (phi !=  Inf) ? phi - delta*(npall-1-n) :
            axis_points[d][npaxis-1] - delta*(npaxis-1-global_id);
      }

      double coord = mesh.get_axis_point(d, id);
      ASSERT_EQ(expval1, coord);
      ASSERT_EQ(expval2, coord);
    }
  }
}


#endif  // WONTON_ENABLE_MPI
