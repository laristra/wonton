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

template<long D>
double refinement_function_2d(const Wonton::Point<D> r) {
  return(1.2 - 1.1*r[0] + 2.4*r[1]);
}

// ============================================================================

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

  // Create a simple 2x2 mesh
  const Wonton::Point<D> lo = {0.0, 0.0};
  const Wonton::Point<D> hi = {1.0, 1.0};
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function_2d, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  ASSERT_EQ(mesh.num_cells(), 46);    // Known from testing

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

}


// ============================================================================

/*TEST(Adaptive_Refinement_Mesh, SmallGrid1D) {

  const int D = 1;

  // Create a single cell mesh
  //                           2^-4    2^-3   2^-2  2^-1 2^0
  const std::vector<double> edges1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  std::vector<double> edges[D];
  for (int d = 0; d < D; ++d) {
    edges[d] = edges1;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  ASSERT_EQ(mesh.axis_num_points(0), edges1.size());

  // Cell coordinates
  for (int d = 0; d < D; ++d) {
    for (int n = 0; n < mesh.axis_num_points(d); ++n) {
      ASSERT_EQ(mesh.axis_point_coordinate(d,n), std::pow(0.5,4-n));
    }
  }

}*/


