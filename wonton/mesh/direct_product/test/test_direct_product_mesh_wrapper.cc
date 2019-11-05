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
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

#include "gtest/gtest.h"

// ============================================================================

namespace direct_product_mesh_wrapper_test {
  template<int D>
  void check_basic_functions(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<std::vector<double>,(std::size_t)D> &axis_points) {

    // Check basic dimensionality
    ASSERT_EQ(mesh_wrapper.space_dimension(), D);

    // Count cells in mesh
    int cell_count = 1;
    for (int d = 0; d < D; ++d) {
      cell_count *= mesh_wrapper.axis_num_cells(d);
    }
    ASSERT_EQ(mesh_wrapper.num_owned_cells(), cell_count);

    // For each axis
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(mesh_wrapper.axis_num_cells(d), (int) axis_points[d].size()-1);
      int n = 0;
      for (auto iter = mesh_wrapper.axis_point_begin(d);
           iter != mesh_wrapper.axis_point_end(d);
           ++iter) {
        ASSERT_EQ(mesh_wrapper.get_axis_point(d,*iter),
            axis_points[d][n]);
        n++;
      }
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void check_cell_bounds(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<int,D>& indices,
      const std::array<std::vector<double>,D> &axis_points) {
    // Get the cell ID
    int id = mesh_wrapper.indices_to_cellid(indices);
    // Get the bounding box
    Wonton::Point<D> plo, phi;
    mesh_wrapper.cell_get_bounds(id, &plo, &phi);
    // Verify the bounding box
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(plo[d], axis_points[d][indices[d]]);
      ASSERT_EQ(phi[d], axis_points[d][indices[d]+1]);
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void check_indices_and_cellids(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<int,D>& indices, const int id) {
    // Ensure indices to cell ID works.
    ASSERT_EQ(mesh_wrapper.indices_to_cellid(indices), id);
    // Ensure cell ID to indices works.
    std::array<int,D> indices2 = mesh_wrapper.cellid_to_indices(id);
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(indices2[d], indices[d]);
    }
    // Ensure cell ID to indices to cell ID works.
    int id2 = mesh_wrapper.indices_to_cellid(
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
      const int d, std::array<int,D> &indices, int &id,
      const Wonton::Direct_Product_Mesh_Wrapper<D> &mesh_wrapper,
      const std::array<std::vector<double>,D> &axis_points) {
    if (d >= 0) {  // Recursive case
      auto num_cells = mesh_wrapper.axis_num_cells(d);
      for (indices[d] = 0; indices[d] < num_cells; ++indices[d]) {
        loop_over_grid<D>(d-1, indices, id, mesh_wrapper, axis_points);
      }
    } else {  // Base case
      // Verify cell bounding boxes
      // -- Since this mesh is just axis-aligned boxes, the bounding boxes will
      //    simply be the cell bounds.
      direct_product_mesh_wrapper_test::check_cell_bounds<D>(
          mesh_wrapper, indices, axis_points);
      // Check IDs vs indices
      direct_product_mesh_wrapper_test::check_indices_and_cellids<D>(
          mesh_wrapper, indices, id);
      // Increment to next cell ID
      id++;
    }
  }
  
  // --------------------------------------------------------------------------

  template<int D>
  void verify_cell_iteration(
      Wonton::Direct_Product_Mesh_Wrapper<D>& wrapper) {
    // Declare an array of zeros for storage
    std::vector<int> found(wrapper.num_owned_cells(), 0);
    // Loop over all cell IDs
    auto ekind = Wonton::Entity_kind::CELL;
    auto etype = Wonton::Entity_type::PARALLEL_OWNED;
    for (auto iter = wrapper.begin(ekind,etype);
        iter != wrapper.end(ekind,etype); ++iter) {
      // Get the indices from the cell ID
      auto id = (*iter);
      auto indices = wrapper.cellid_to_indices(id);
      // Convert the indices to an arbitrary 1D index for storage array
      std::size_t n = indices[0];
      for (int d = 1; d < D; ++d) {
        n *= wrapper.axis_num_cells(d-1);
        n += indices[d];
      }
      // Ensure that the 1D index is in the valid range
      ASSERT_GE(n, 0);
      ASSERT_LT(n, wrapper.num_owned_cells());
      // Mark that value as having been found
      found[n] = 1;
    }
    // Verify that all elements were found
    auto num_found = std::accumulate(found.begin(), found.end(), 0);
    ASSERT_EQ(num_found, wrapper.num_owned_cells());
  }

  // --------------------------------------------------------------------------

  template<int D>
  void run_all_tests(const std::vector<double> &axis_points1) {
    // Build mesh and wrapper
    std::array<std::vector<double>,D> axis_points;
    for (int d = 0; d < D; ++d) {
      axis_points[d] = axis_points1;
    }
    Wonton::Direct_Product_Mesh<D> mesh(axis_points);
    Wonton::Direct_Product_Mesh_Wrapper<D> mesh_wrapper(mesh);

    // Run basic tests
    direct_product_mesh_wrapper_test::check_basic_functions(
        mesh_wrapper, axis_points);

    // Run looping tests
    int id = 0;
    std::array<int,D> indices;
    direct_product_mesh_wrapper_test::loop_over_grid<D>(
        D-1, indices, id, mesh_wrapper, axis_points);

    // Verify cell iteration
    verify_cell_iteration(mesh_wrapper);
  }

} // namespace direct_product_mesh_wrapper_test 


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid1D) {
  const std::vector<double> axis_points1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<1>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid2D) {
  const std::vector<double> axis_points1 = {0.0, 0.25, 0.75, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<2>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, OneCell3D) {
  const std::vector<double> axis_points1 = {0.0, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<3>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmalGrid10D) {
  const std::vector<double> axis_points1 = {0.0, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<10>(axis_points1);
}


// ============================================================================

