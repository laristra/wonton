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
TEST(Direct_Product_Mesh, SingleCell3D) {

  const int D = 3;

  // Create a single cell mesh
  const std::vector<double> edges1 = {0.0, 1.0};
  std::vector<double> edges[D];
  for (int d = 0; d < D; ++d) {
    edges[d] = edges1;
  }
  Wonton::Direct_Product_Mesh<D> mesh(edges);

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

}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid2D) {

  const int D = 2;

  // Create a simple 2x2 mesh
  const std::vector<double> edges1 = {0.0, 0.5, 1.0};
  std::vector<double> edges[D];
  for (int d = 0; d < D; ++d) {
    edges[d] = edges1;
  }
  Wonton::Direct_Product_Mesh<D> mesh(edges);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  for (int d = 0; d < D; ++d) {
    ASSERT_EQ(mesh.axis_num_points(d), edges1.size());
  }

  // Cell coordinates
  for (int d = 0; d < D; ++d) {
    for (int n = 0; n < mesh.axis_num_points(d); ++n) {
      ASSERT_EQ(mesh.axis_point_coordinate(d,n), 0.0 + 0.5*n);
    }
  }

}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid1D) {

  const int D = 1;

  // Create a single cell mesh
  //                           2^-4    2^-3   2^-2  2^-1 2^0
  const std::vector<double> edges1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  std::vector<double> edges[D];
  for (int d = 0; d < D; ++d) {
    edges[d] = edges1;
  }
  Wonton::Direct_Product_Mesh<D> mesh(edges);

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

}


