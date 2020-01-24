/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#include <array>
#include <iostream>
#include <iterator>
#include <numeric>
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
    ASSERT_TRUE(mesh_wrapper.num_owned_cells() >= 0);
    ASSERT_TRUE(mesh_wrapper.num_ghost_cells() >= 0);
    ASSERT_EQ(cell_count, mesh_wrapper.num_owned_cells()+mesh_wrapper.num_ghost_cells());

    // For each axis check the points
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(mesh_wrapper.axis_num_cells(d), (int) axis_points[d].size()-1);

      // We have to ignore the first and last point as they are ghosts
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
    if (mesh_wrapper.distributed()) {
      // axis_points[0] is the ghost point (index = -1) while
      // axis_points[1] is the first owned point (index = 0)
      for (int d = 0; d < D; ++d) {
        ASSERT_EQ(plo[d], axis_points[d][indices[d]+1]);
        ASSERT_EQ(phi[d], axis_points[d][indices[d]+2]);
      }
    } else {  // serial case
      // axis_points[0] is first owned point (index = 0);
      for (int d = 0; d < D; ++d) {
        ASSERT_EQ(plo[d], axis_points[d][indices[d]]);
        ASSERT_EQ(phi[d], axis_points[d][indices[d]+1]);
      }
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
      auto num_cells_all = mesh_wrapper.axis_num_cells(d);
      if (mesh_wrapper.distributed()) {
        int start_index = -1;
        for (indices[d] = start_index; indices[d] < num_cells_all-1;
             ++indices[d]) {
          loop_over_grid<D>(d-1, indices, id, mesh_wrapper, axis_points);
        }
      } else {
        int start_index = 0;
        for (indices[d] = start_index; indices[d] < num_cells_all;
             ++indices[d]) {
          loop_over_grid<D>(d-1, indices, id, mesh_wrapper, axis_points);
        }
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
    int num_cells_all = wrapper.num_owned_cells()+wrapper.num_ghost_cells();
    std::vector<int> found(num_cells_all, 0);
    // Loop over all cell IDs
    auto ekind = Wonton::Entity_kind::CELL;
    auto etype = Wonton::Entity_type::ALL;
    for (auto iter = wrapper.begin(ekind,etype);
        iter != wrapper.end(ekind,etype); ++iter) {
      // Get the indices from the cell ID
      auto id = (*iter);
      auto indices = wrapper.cellid_to_indices(id);
      // Convert the indices to an arbitrary 1D index for storage array
      std::size_t n = 0;
      for (int d = D-1; d > 0; --d) {
        n += (wrapper.distributed() ? indices[d]+1 : indices[d]);
        n *= wrapper.axis_num_cells(d-1);
      }
      n += wrapper.distributed() ? indices[0]+1 : indices[0];
      // Ensure that the 1D index is in the valid range
      ASSERT_GE(n, 0);
      ASSERT_LT(n, num_cells_all);
      // Mark that value as having been found
      found[n] += 1;
    }
    // Verify that all elements were found
    // -- Did we find the right number?
    auto num_found = std::accumulate(found.begin(), found.end(), 0);
    ASSERT_EQ(num_found, num_cells_all);
    // -- Was every element found?
    int minimum = *std::min_element(found.begin(), found.end());
    ASSERT_EQ(minimum, 1);
    // -- Was every element found only once?
    int maximum = *std::max_element(found.begin(), found.end());
    ASSERT_EQ(maximum, 1);
  }

  // --------------------------------------------------------------------------

  template<int D>
  void run_all_tests(const std::vector<double> &axis_points_in) {
    // Build mesh and wrapper
    std::array<std::vector<double>,D> axis_points_global;
    for (int d = 0; d < D; ++d) {
      axis_points_global[d] = axis_points_in;
    }


    Wonton::Executor_type *executor;
#ifdef WONTON_ENABLE_MPI
    Wonton::MPIExecutor_type parallel_exec(MPI_COMM_WORLD);
    executor = &parallel_exec;
#else
    Wonton::SerialExecutor_type serial_exec;
    executor = &serial_exec;
#endif
    
    Wonton::Direct_Product_Mesh<D> mesh(axis_points_global, executor);
    Wonton::Direct_Product_Mesh_Wrapper<D> mesh_wrapper(mesh);

    // collect the points on this partition (all points bounding OWNED
    // cells, which means all points except the OUTER two in any
    // direction)
    std::array<std::vector<double>,D> axis_points;
    for (int d = 0; d < D; d++) {
      int n = mesh.axis_num_points(d, Wonton::ALL);
      axis_points[d].resize(n);
      if (mesh.distributed())
        for (int i = -1; i < n-1; i++)
          axis_points[d][i+1] = mesh.get_axis_point(d, i);
      else
        for (int i = 0; i < n; i++)
          axis_points[d][i] = mesh.get_axis_point(d, i);
    }
    
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

// TEST(Direct_Product_Mesh_Wrapper, SmallGrid1D) {
//   const std::vector<double> axis_points1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
//   direct_product_mesh_wrapper_test::run_all_tests<1>(axis_points1);
// }


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

TEST(Direct_Product_Mesh_Wrapper, SmallGrid10D) {
  const std::vector<double> axis_points1 = {0.0, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<10>(axis_points1);
}


// ============================================================================

