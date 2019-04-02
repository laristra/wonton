/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

// ============================================================================
/*!
  @file test_adaptive_refinement_mesh.cc
  @brief Tests for the adaptive_refinement mesh definition in
         adaptive_refinement_mesh.h
*/

// ============================================================================

double refinement_function_1d(const Wonton::Point<1> r) {
  return(2.0 + 3.0*r[0]);
}

// ============================================================================

double refinement_function_2d(const Wonton::Point<2> r) {
  return(1.2 - 1.1*r[0] + 2.4*r[1]);
}

// ============================================================================

// TODO: 3D test
/*TEST(Adaptive_Refinement_Mesh, SingleCell3D) {

  const int D = 3;

  // Create a single cell mesh
  const std::vector<double> edges1 = {0.0, 1.0};
  std::vector<double> edges[D];
  for (int d = 0; d < D; ++d) {
    edges[d] = edges1;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  for (int d = 0; d < D; ++d) {
    ASSERT_EQ(mesh.axis_num_points(d), edges1.size());
  }

  // Cell coordinates
  for (int d = 0; d < D; ++d) {
    for (int n = 0; n < mesh.axis_num_points(d); ++n) {
      ASSERT_EQ(mesh.axis_point_coordinate(d,n), 0.0 + 1.0*n);
    }
  }

}*/


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test2D) {

  const int D = 2;

  // Create a mesh
  const Wonton::Point<D> lo = {0.0, 0.0};
  const Wonton::Point<D> hi = {1.0, 1.0};
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function_2d, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 46);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(4);
  ASSERT_EQ(bounds[0][0], 0.5);
  ASSERT_EQ(bounds[0][1], 1.0);
  ASSERT_EQ(bounds[1][0], 0.0);
  ASSERT_EQ(bounds[1][1], 0.5);
  bounds = mesh.cell_get_bounds(20);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 1.0/16.0);
  ASSERT_EQ(bounds[1][0], 1.0-1.0/16.0);
  ASSERT_EQ(bounds[1][1], 1.0);
  bounds = mesh.cell_get_bounds(35);
  ASSERT_EQ(bounds[0][0], 0.5);
  ASSERT_EQ(bounds[0][1], 5.0/8.0);
  ASSERT_EQ(bounds[1][0], 5.0/8.0);
  ASSERT_EQ(bounds[1][1], 3.0/4.0);

}


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test1D) {

  const int D = 1;

  // Create a mesh
  const Wonton::Point<D> lo = {0.0};
  const Wonton::Point<D> hi = {1.0};
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function_1d, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 18);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(0);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 0.125);
  bounds = mesh.cell_get_bounds(6);
  ASSERT_EQ(bounds[0][0], 9.0/16.0);
  ASSERT_EQ(bounds[0][1], 10.0/16.0);
  bounds = mesh.cell_get_bounds(14);
  ASSERT_EQ(bounds[0][0], 1.0 - 4.0/32.0);
  ASSERT_EQ(bounds[0][1], 1.0 - 3.0/32.0);

}


