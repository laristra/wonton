/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#include <array>
#include <iostream>
#include <iterator>
#include <vector>

#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"
#include "wonton/support/wonton.h"
#include "wonton/support/CellID.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

#include "gtest/gtest.h"

// ============================================================================

namespace direct_product_mesh_wrapper_test {
  template<int D>
  void check_basic_functions(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<std::vector<double>,(std::size_t)D> &edges) {

    // Check basic dimensionality
    ASSERT_EQ(mesh_wrapper.space_dimension(), D);

    // Count cells in mesh
    int cell_count = 1;
    for (int d = 0; d < D; ++d) {
      cell_count *= mesh_wrapper.axis_num_cells(d);
    }
    ASSERT_EQ(mesh_wrapper.total_num_cells(), cell_count);

    // For each axis
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(mesh_wrapper.axis_num_cells(d), edges[d].size()-1);
      int n = 0;
      for (auto iter = mesh_wrapper.axis_point_begin(d);
           iter != mesh_wrapper.axis_point_end(d);
           ++iter) {
        ASSERT_EQ(mesh_wrapper.axis_point_coordinate(d,*iter), edges[d][n]);
        n++;
      }
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void check_cell_bounds(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<int,D>& indices,
      const std::array<std::vector<double>,D> &edges) {
    // Get the cell ID
    Wonton::CellID id = mesh_wrapper.indices_to_cellid(indices);
    // Get the bounding box
    Wonton::Point<D> plo, phi;
    mesh_wrapper.cell_get_bounds(id, &plo, &phi);
    // Verify the bounding box
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(plo[d], edges[d][indices[d]]);
      ASSERT_EQ(phi[d], edges[d][indices[d]+1]);
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void check_indices_and_cellids(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<int,D>& indices, const Wonton::CellID id) {
    // Ensure indices to cell ID works.
    ASSERT_EQ(mesh_wrapper.indices_to_cellid(indices), id);
    // Ensure cell ID to indices works.
    std::array<int,D> indices2 = mesh_wrapper.cellid_to_indices(id);
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(indices2[d], indices[d]);
    }
    // Ensure cell ID to indices to cell ID works.
    Wonton::CellID id2 = mesh_wrapper.indices_to_cellid(
        mesh_wrapper.cellid_to_indices(id));
    ASSERT_EQ(id2, id);
    // Ensure indices to cell ID to indices works.
    indices2 = mesh_wrapper.cellid_to_indices(
        mesh_wrapper.indices_to_cellid(indices));
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(indices2[d], indices[d]);
    }
  }
  
  // --------------------------------------------------------------------------

  template<int D>
  void loop_over_grid(
      const int d, std::array<int,D> &indices, Wonton::CellID &id,
      const Wonton::Direct_Product_Mesh_Wrapper<D> &mesh_wrapper,
      const std::array<std::vector<double>,D> &edges) {
    if (d >= 0) {  // Recursive case
      auto num_cells = mesh_wrapper.axis_num_cells(d);
      for (indices[d] = 0; indices[d] < num_cells; ++indices[d]) {
        loop_over_grid<D>(d-1, indices, id, mesh_wrapper, edges);
      }
    } else {  // Base case
      // Verify cell bounding boxes
      // -- Since this mesh is just axis-aligned boxes, the bounding boxes will
      //    simply be the cell bounds.
      direct_product_mesh_wrapper_test::check_cell_bounds<D>(
          mesh_wrapper, indices, edges);
      // Check IDs vs indices
      direct_product_mesh_wrapper_test::check_indices_and_cellids<D>(
          mesh_wrapper, indices, id);
      // Increment to next cell ID
      id++;
    }
  }
  
  // --------------------------------------------------------------------------

  template<int D>
  void run_all_tests(const std::vector<double> &edges1) {
    // Build mesh and wrapper
    std::array<std::vector<double>,D> edges;
    for (int d = 0; d < D; ++d) {
      edges[d] = edges1;
    }
    Wonton::Direct_Product_Mesh<D> mesh(edges);
    Wonton::Direct_Product_Mesh_Wrapper<D> mesh_wrapper(mesh);

    // Run basic tests
    direct_product_mesh_wrapper_test::check_basic_functions(
        mesh_wrapper, edges);

    // Run looping tests
    Wonton::CellID id = 0;
    std::array<int,D> indices;
    direct_product_mesh_wrapper_test::loop_over_grid<D>(
        D-1, indices, id, mesh_wrapper, edges);
  }

} // namespace direct_product_mesh_wrapper_test 


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid1D) {
  const std::vector<double> edges1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<1>(edges1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid2D) {
  const std::vector<double> edges1 = {0.0, 0.25, 0.75, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<2>(edges1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, OneCell3D) {
  const std::vector<double> edges1 = {0.0, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<3>(edges1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmalGrid10D) {
  const std::vector<double> edges1 = {0.0, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<10>(edges1);
}


// ============================================================================

