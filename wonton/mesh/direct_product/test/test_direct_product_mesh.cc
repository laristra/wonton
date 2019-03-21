/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/direct_product/direct_product_mesh.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

// ============================================================================
/*!
  @file test_direct_product_mesh.cc
  @brief Tests for the direct_product mesh definition in direct_product_mesh.h
*/
TEST(Direct_Product_Mesh3D, SingleCell) {

  // Create a single cell mesh
  std::vector<double> edges = {0.0, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges, edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 3);

  // Cell counts
  for (int dim = 0; dim < 3; ++dim) {
    ASSERT_EQ(mesh.axis_num_points(dim), 2);
  }

  // Cell coordinates
  for (int dim = 0; dim < 3; ++dim) {
    for (int n = 0; n < mesh.axis_num_points(dim); ++n) {
      ASSERT_EQ(mesh.axis_point_coordinate(dim, n), 0.0 + 1.0*n);
    }
  }

}


// ============================================================================

TEST(Direct_Product_Mesh2D, SmallGrid) {

  // Create a simple 2x2 mesh
  std::vector<double> edges = {0.0, 0.5, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 2);

  // Cell counts
  for (int dim = 0; dim < 2; ++dim) {
    ASSERT_EQ(mesh.axis_num_points(dim), 3);
  }

  // Cell coordinates
  for (int dim = 0; dim < 2; ++dim) {
    for (int n = 0; n < mesh.axis_num_points(dim); ++n) {
      ASSERT_EQ(mesh.axis_point_coordinate(dim, n), 0.0 + 0.5*n);
    }
  }

}


// ============================================================================

TEST(Direct_Product_Mesh1D, SmallGrid) {

  // Create a single cell mesh
  //                           2^-4    2^-3   2^-2  2^-1 2^0
  std::vector<double> edges = {0.0625, 0.125, 0.25, 0.5, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 1);

  // Cell counts
  ASSERT_EQ(mesh.axis_num_points(0), 5);

  // Cell coordinates
  for (int n = 0; n < mesh.axis_num_points(0); ++n) {
    ASSERT_EQ(mesh.axis_point_coordinate(0, n), std::pow(0.5,4-n));
  }

}


