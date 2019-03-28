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

namespace direct_product_mesh_wrapper_test {
  template<long dim>
  void check_basic_functions(
      const Wonton::Direct_Product_Mesh_Wrapper& mesh_wrapper,
      const std::vector<double>& edges) {

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
  }

  // --------------------------------------------------------------------------

  template<long dim>
  void check_cell_bounds(
      const Wonton::Direct_Product_Mesh_Wrapper& mesh_wrapper,
      const Wonton::IntPoint<dim>& indices,
      const std::vector<double>& edges) {
    // Get the cell ID
    Wonton::CellID id = mesh_wrapper.indices_to_cellid<dim>(indices);
    // Get the bounding box
    Wonton::Point<dim> plo, phi;
    mesh_wrapper.cell_get_bounds(id, &plo, &phi);
    // Verify the bounding box (assumes edges are same for all axes)
    for (int d = 0; d < dim; ++d) {
      ASSERT_EQ(plo[d], edges[indices[d]]);
      ASSERT_EQ(phi[d], edges[indices[d]+1]);
    }
  }

  // --------------------------------------------------------------------------

  template<long dim>
  void check_indices_and_cellids(
      const Wonton::Direct_Product_Mesh_Wrapper& mesh_wrapper,
      const Wonton::IntPoint<dim>& indices, const Wonton::CellID id) {
    // Ensure indices to cell ID works.
    ASSERT_EQ(mesh_wrapper.indices_to_cellid<dim>(indices), id);
    // Ensure cell ID to indices works.
    Wonton::IntPoint<dim> indices2 = mesh_wrapper.cellid_to_indices<dim>(id);
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
  }

} // namespace direct_product_mesh_wrapper_test 

// ============================================================================

TEST(Direct_Product_Mesh, OneCell3D) {

  const int dim = 3;

  // Build a mesh and wrapper
  std::vector<double> edges = {0.0, 1.0};
  Wonton::Direct_Product_Mesh mesh(edges, edges, edges);
  Wonton::Direct_Product_Mesh_Wrapper mesh_wrapper(mesh);

  // Run basic tests
  direct_product_mesh_wrapper_test::check_basic_functions<dim>(
      mesh_wrapper, edges);

  Wonton::CellID id = 0;
  for (int k = 0; k < mesh_wrapper.axis_num_cells(2); ++k) {
    for (int j = 0; j < mesh_wrapper.axis_num_cells(1); ++j) {
      for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
        const Wonton::IntPoint<dim> indices = {i, j, k};
        // Verify cell bounding boxes
        // -- Since this mesh is just axis-aligned boxes, the bounding boxes
        //    will simply be the cell bounds.
        direct_product_mesh_wrapper_test::check_cell_bounds<dim>(
            mesh_wrapper, indices, edges);
        // Check IDs vs indices
        direct_product_mesh_wrapper_test::check_indices_and_cellids<dim>(
            mesh_wrapper, indices, id);
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

  // Run basic tests
  direct_product_mesh_wrapper_test::check_basic_functions<dim>(
      mesh_wrapper, edges);

  Wonton::CellID id = 0;
  for (int j = 0; j < mesh_wrapper.axis_num_cells(1); ++j) {
    for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
      const Wonton::IntPoint<dim> indices = {i, j};
      // Verify cell bounding boxes
      // -- Since this mesh is just axis-aligned boxes, the bounding boxes will
      //    simply be the cell bounds.
      direct_product_mesh_wrapper_test::check_cell_bounds<dim>(
          mesh_wrapper, indices, edges);
      // Check IDs vs indices
      direct_product_mesh_wrapper_test::check_indices_and_cellids<dim>(
          mesh_wrapper, indices, id);
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

  // Run basic tests
  direct_product_mesh_wrapper_test::check_basic_functions<dim>(
      mesh_wrapper, edges);

  Wonton::CellID id = 0;
  for (int i = 0; i < mesh_wrapper.axis_num_cells(0); ++i) {
    const Wonton::IntPoint<dim> indices = {i};
    // Verify cell bounding boxes
    // -- Since this mesh is just axis-aligned boxes, the bounding boxes will
    //    simply be the cell bounds.
    direct_product_mesh_wrapper_test::check_cell_bounds<dim>(
        mesh_wrapper, indices, edges);
    // Check IDs vs indices
    direct_product_mesh_wrapper_test::check_indices_and_cellids<dim>(
        mesh_wrapper, indices, id);
    // Increment to next cell ID
    id++;
  }

}


// ============================================================================

