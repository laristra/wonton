/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#include <iostream>
#include <vector>
#include <iterator>

#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

#include "gtest/gtest.h"

// ============================================================================

TEST(Direct_Product_Mesh, OneCell3D) {

  const int dim = 3;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges, edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Build a useful structure
  std::vector<double> all_edges[dim];
  for (int d = 0; d < dim; ++d) {
    all_edges[dim] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[dim][n]);
      n++;
    }
  }
}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid2D) {

  const int dim = 2;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0, 0.25, 0.75, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges, edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Build a useful structure
  std::vector<double> all_edges[dim];
  for (int d = 0; d < dim; ++d) {
    all_edges[dim] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[dim][n]);
      n++;
    }
  }
}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid1D) {

  const int dim = 1;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0625, 0.125, 0.25, 0.5, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges, edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Build a useful structure
  std::vector<double> all_edges[dim];
  for (int d = 0; d < dim; ++d) {
    all_edges[dim] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[dim][n]);
      n++;
    }
  }
}


// ============================================================================

