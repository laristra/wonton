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
#include "wonton/support/CellID.h"
#include "wonton/support/IntPoint.h"
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
    all_edges[d] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // Count cells in mesh
  int cell_count = 1;
  for (int d = 0; d < dim; ++d) {
    cell_count *= mesh_wrapper.axis_num_cells(d);
  }
  ASSERT_EQ(mesh_wrapper.total_num_cells(), cell_count);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[d][n]);
      n++;
    }
  }

  // Check IDs vs indices
  Wonton::CellID id = 0;
  for (int k = 0; k < mesh_wrapper.axis_num_cells(2); ++k) {
    for (int j = 0; j < mesh_wrapper.axis_num_cells(1); ++j) {
      for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
        const Wonton::IntPoint<dim> indices = {i, j, k};
        // Ensure indices to cell ID works.
        ASSERT_EQ(mesh_wrapper.indices_to_cellid<dim>(indices), id);
        Wonton::IntPoint<dim> indices2 = 
          mesh_wrapper.cellid_to_indices<dim>(id);
        // Ensure cell ID to indices works.
        for (int d = 0; d < dim; ++d) {
          ASSERT_EQ(indices2[d], indices[d]);
        }
        // Ensure cell ID to indices to cell ID works.
        Wonton::CellID id2 = mesh_wrapper.indices_to_cellid<dim>(
            mesh_wrapper.cellid_to_indices<dim>(id));
        ASSERT_EQ(id2, id);
        // Ensure indices to cell ID to indices works.
        indices2 = mesh_wrapper.cellid_to_indices<dim>(
            mesh_wrapper.indices_to_cellid<dim>(indices));
        for (int d = 0; d < dim; ++d) {
          ASSERT_EQ(indices2[d], indices[d]);
        }
        // Increment to next cell ID
        id++;
      }
    }
  }

}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid2D) {

  const int dim = 2;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0, 0.25, 0.75, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Build a useful structure
  std::vector<double> all_edges[dim];
  for (int d = 0; d < dim; ++d) {
    all_edges[d] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // Count cells in mesh
  int cell_count = 1;
  for (int d = 0; d < dim; ++d) {
    cell_count *= mesh_wrapper.axis_num_cells(d);
  }
  ASSERT_EQ(mesh_wrapper.total_num_cells(), cell_count);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[d][n]);
      n++;
    }
  }

  // Check IDs vs indices
  Wonton::CellID id = 0;
  for (int j = 0; j < mesh_wrapper.axis_num_cells(1); ++j) {
    for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
      const Wonton::IntPoint<dim> indices = {i, j};
      // Ensure indices to cell ID works.
      ASSERT_EQ(mesh_wrapper.indices_to_cellid<dim>(indices), id);
      Wonton::IntPoint<dim> indices2 = mesh_wrapper.cellid_to_indices<dim>(id);
      // Ensure cell ID to indices works.
      for (int d = 0; d < dim; ++d) {
        ASSERT_EQ(indices2[d], indices[d]);
      }
      // Ensure cell ID to indices to cell ID works.
      Wonton::CellID id2 = mesh_wrapper.indices_to_cellid<dim>(
          mesh_wrapper.cellid_to_indices<dim>(id));
      ASSERT_EQ(id2, id);
      // Ensure indices to cell ID to indices works.
      indices2 = mesh_wrapper.cellid_to_indices<dim>(
        mesh_wrapper.indices_to_cellid<dim>(indices));
      for (int d = 0; d < dim; ++d) {
        ASSERT_EQ(indices2[d], indices[d]);
      }
      // Increment to next cell ID
      id++;
    }
  }

}


// ============================================================================

TEST(Direct_Product_Mesh, SmallGrid1D) {

  const int dim = 1;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0625, 0.125, 0.25, 0.5, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Build a useful structure
  std::vector<double> all_edges[dim];
  for (int d = 0; d < dim; ++d) {
    all_edges[d] = edges;
  }

  // Check basic dimensionality
  ASSERT_EQ(mesh_wrapper.space_dimension(), dim);

  // Count cells in mesh
  int cell_count = 1;
  for (int d = 0; d < dim; ++d) {
    cell_count *= mesh_wrapper.axis_num_cells(d);
  }
  ASSERT_EQ(mesh_wrapper.total_num_cells(), cell_count);

  // For each axis
  for (int d = 0; d < dim; ++d) {
    ASSERT_EQ(mesh_wrapper.axis_num_cells(d), all_edges[d].size()-1);
    int n = 0;
    for (auto iter = mesh_wrapper.axis_point_begin(d);
         iter != mesh_wrapper.axis_point_end(d);
         ++iter) {
      ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), all_edges[d][n]);
      n++;
    }
  }

  // Check IDs vs indices
  Wonton::CellID id = 0;
  for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
    const Wonton::IntPoint<dim> indices = {i};
    // Ensure indices to cell ID works.
    ASSERT_EQ(mesh_wrapper.indices_to_cellid<dim>(indices), id);
    Wonton::IntPoint<dim> indices2 = mesh_wrapper.cellid_to_indices<dim>(id);
    // Ensure cell ID to indices works.
    for (int d = 0; d < dim; ++d) {
      ASSERT_EQ(indices2[d], indices[d]);
    }
    // Ensure cell ID to indices to cell ID works.
    Wonton::CellID id2 = mesh_wrapper.indices_to_cellid<dim>(
        mesh_wrapper.cellid_to_indices<dim>(id));
    ASSERT_EQ(id2, id);
    // Ensure indices to cell ID to indices works.
    indices2 = mesh_wrapper.cellid_to_indices<dim>(
      mesh_wrapper.indices_to_cellid<dim>(indices));
    for (int d = 0; d < dim; ++d) {
      ASSERT_EQ(indices2[d], indices[d]);
    }
    // Increment to next cell ID
    id++;
  }

}


// ============================================================================

