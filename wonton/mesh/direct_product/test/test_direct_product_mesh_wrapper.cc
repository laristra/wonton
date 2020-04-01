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

    // Is the mesh distributed?
    bool is_distributed;
#ifdef WONTON_ENABLE_MPI
    is_distributed = true;
#else
    is_distributed = false;
#endif
    ASSERT_EQ(mesh_wrapper.distributed(), is_distributed);

    // Check mesh bounds
    Point<D> plo, phi;
    mesh_wrapper.get_global_bound(&plo, &phi);
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(plo[d], axis_points[d].front());
      ASSERT_EQ(phi[d], axis_points[d].back());
    }

    // Count cells in mesh
    int cell_count = 1;
    for (int d = 0; d < D; ++d) {
      cell_count *= mesh_wrapper.num_axis_cells(d);
    }
    ASSERT_GE(mesh_wrapper.num_owned_cells(), 0);
    ASSERT_GE(mesh_wrapper.num_ghost_cells(), 0);
    ASSERT_EQ(cell_count,
              mesh_wrapper.num_owned_cells()+mesh_wrapper.num_ghost_cells());

    // For each axis check the points
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(mesh_wrapper.num_axis_cells(d), (int) axis_points[d].size()-1);

      // We have to ignore the first and last point as they are ghosts
      int n = 0;
      for (auto iter = mesh_wrapper.axis_point_begin(d);
           iter != mesh_wrapper.axis_point_end(d);
           ++iter) {
        ASSERT_EQ(mesh_wrapper.get_axis_point(d,*iter), axis_points[d][n]);
        n++;
      }
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void check_cell_geometry(
      const Wonton::Direct_Product_Mesh_Wrapper<D>& mesh_wrapper,
      const std::array<int,D>& indices,
      const std::array<std::vector<double>,D> &axis_points) {
    int num_ghost_layers = mesh_wrapper.num_ghost_layers();
    // Get the cell ID
    int id = mesh_wrapper.indices_to_cellid(indices);
    // Get the bounding box
    Wonton::Point<D> plo, phi;
    mesh_wrapper.cell_get_bounds(id, &plo, &phi);

    // Verify the bounding box
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(plo[d], axis_points[d][indices[d]+num_ghost_layers]);
      ASSERT_EQ(phi[d], axis_points[d][indices[d]+num_ghost_layers+1]);
    }

    // Verify cell coordinate count
    std::vector<Wonton::Point<D>> ccoord;
    mesh_wrapper.cell_get_coordinates(id, &ccoord);
    ASSERT_EQ(pow(2,D), ccoord.size());
      
    if (mesh_wrapper.cell_get_type(id) == Wonton::PARALLEL_OWNED) {
      // Verify the volume
      double expvolume = 1.0;
      for (int d = 0; d < D; d++) {
        double lo = axis_points[d][indices[d]+num_ghost_layers];
        double hi = axis_points[d][indices[d]+num_ghost_layers+1];
        expvolume *= (hi-lo);
      }
      double volume = mesh_wrapper.cell_volume(id);
      ASSERT_EQ(expvolume, volume);
      
      // Verify the cell centroid
      Wonton::Point<D> expcen;
      for (int d = 0; d < D; d++) expcen[d] = 0.0;
      for (int i = 0; i < ccoord.size(); i++)
        expcen += ccoord[i];
      expcen /= ccoord.size();
      Wonton::Point<D> ccen;
      mesh_wrapper.cell_centroid(id, &ccen);
      for (int d = 0; d < D; d++)
        ASSERT_NEAR(expcen[d], ccen[d], 1.0e-14);
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
      auto num_cells_all = mesh_wrapper.num_axis_cells(d);
      int num_ghost_layers = mesh_wrapper.num_ghost_layers();
      for (int i = 0; i < num_cells_all; i++) {
        indices[d] = i - num_ghost_layers;
        loop_over_grid<D>(d-1, indices, id, mesh_wrapper, axis_points);
      }
    } else {  // Base case
      // Verify cell bounding boxes
      // -- Since this mesh is just axis-aligned boxes, the bounding boxes will
      //    simply be the cell bounds.
      direct_product_mesh_wrapper_test::check_cell_geometry<D>(
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
    Wonton::Entity_kind ekind = Wonton::CELL;
    Wonton::Entity_type etype = Wonton::ALL;
    for (auto iter = wrapper.begin(ekind,etype);
        iter != wrapper.end(ekind,etype); ++iter) {
      // Get the indices from the cell ID
      auto id = (*iter);
      auto indices = wrapper.cellid_to_indices(id);
      // Convert the indices to an arbitrary 1D index for storage array
      std::size_t n = 0;
      for (int d = D-1; d > 0; --d) {
        n += indices[d]+wrapper.num_ghost_layers();
        n *= wrapper.num_axis_cells(d-1);
      }
      n += indices[0]+wrapper.num_ghost_layers();
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
  void verify_cell_adjacencies(
      Wonton::Direct_Product_Mesh_Wrapper<D>& wrapper) {    
    // Loop over all cell IDs
    Wonton::Entity_kind ekind = Wonton::CELL;
    Wonton::Entity_type etype = Wonton::ALL;
    for (auto iter = wrapper.begin(ekind,etype);
        iter != wrapper.end(ekind,etype); ++iter) {
      // Get the indices from the cell ID
      auto id = (*iter);
      auto indices = wrapper.cellid_to_indices(id);

      std::vector<int> adjcells;
      wrapper.cell_get_node_adj_cells(id, etype, &adjcells);

      // Make sure we got the right number?
      int num_ghost_layers = wrapper.mesh().num_ghost_layers();
      std::array<int, D> nindices;
      for (int d = 0; d < D; d++) {
        int num_points_all = wrapper.mesh().num_axis_points(d, Wonton::ALL);
        int lo = (indices[d]-1 >= -num_ghost_layers) ?
            indices[d]-1 : indices[d];
        int hi = (indices[d]+1 < num_points_all-num_ghost_layers-1) ?
            indices[d]+1 : indices[d];
        nindices[d] = hi-lo+1;
      }
      int expnum = 1;
      for (int d = 0; d < D; d++)
        expnum *= nindices[d];
      ASSERT_EQ(expnum, adjcells.size());
      
      // Make sure the indices of the cells are no more than one off
      for (auto adjid : adjcells) {
        auto adjindices = wrapper.cellid_to_indices(adjid);
        for (int d = 0; d < D; d++)
          ASSERT_LE(fabs(indices[d]-adjindices[d]), 1);
      }

      bool on_ext_boundary = false;
      for (int d = 0; d < D; d++)
        if (wrapper.mesh().point_on_external_boundary(d, indices[d]) ||
            wrapper.mesh().point_on_external_boundary(d, indices[d]+1)) {
          on_ext_boundary = true;
          break;
        }
      ASSERT_EQ(on_ext_boundary, wrapper.on_exterior_boundary(Wonton::CELL, id));
    }
  }

  // --------------------------------------------------------------------------

  template<int D>
  void verify_node_adjacencies(
      Wonton::Direct_Product_Mesh_Wrapper<D>& wrapper) {    
    // Loop over all node IDs
    Wonton::Entity_kind ekind = Wonton::NODE;
    Wonton::Entity_type etype = Wonton::ALL;
    for (auto iter = wrapper.begin(ekind,etype);
        iter != wrapper.end(ekind,etype); ++iter) {
      // Get the indices from the node ID
      auto id = (*iter);
      auto indices = wrapper.nodeid_to_indices(id);

      std::vector<int> adjnodes;
      wrapper.node_get_cell_adj_nodes(id, etype, &adjnodes);

      // Make sure we got the right number?
      int num_ghost_layers = wrapper.mesh().num_ghost_layers();
      std::array<int, D> nindices;
      for (int d = 0; d < D; d++) {
        int num_points_all = wrapper.mesh().num_axis_points(d, Wonton::ALL);
        int lo = (indices[d]-1 >= -num_ghost_layers) ?
            indices[d]-1 : indices[d];
        int hi = (indices[d]+1 < num_points_all-num_ghost_layers) ?
            indices[d]+1 : indices[d];
        nindices[d] = hi-lo+1;
      }
      int expnum = 1;
      for (int d = 0; d < D; d++)
        expnum *= nindices[d];
      ASSERT_EQ(expnum, adjnodes.size());
      
      // Make sure the indices of the cells are no more than one off
      for (auto adjid : adjnodes) {
        auto adjindices = wrapper.nodeid_to_indices(adjid);
        for (int d = 0; d < D; d++)
          ASSERT_LE(fabs(indices[d]-adjindices[d]), 1);
      }

      bool on_ext_boundary = false;
      for (int d = 0; d < D; d++)
        if (wrapper.mesh().point_on_external_boundary(d, indices[d])) {
          on_ext_boundary = true;
          break;
        }
      ASSERT_EQ(on_ext_boundary, wrapper.on_exterior_boundary(Wonton::NODE, id));      
    }
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
    bool distributed = false;
#ifdef WONTON_ENABLE_MPI
    if (D > 3) return;  // cannot partition in anything greater than 3D
    Wonton::MPIExecutor_type parallel_exec(MPI_COMM_WORLD);
    executor = &parallel_exec;
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs > 1) distributed = true;
#else
    Wonton::SerialExecutor_type serial_exec;
    executor = &serial_exec;
#endif
    
    int maxNG = distributed ? 2 : 0;
    for (int NG = 0; NG <= maxNG; NG++) {
      Wonton::Direct_Product_Mesh<D> mesh(axis_points_global, executor, NG);
      Wonton::Direct_Product_Mesh_Wrapper<D> mesh_wrapper(mesh);

      // collect the points on this partition (all points bounding OWNED
      // cells, which means all points except the OUTER 2*NGhostLayers
      // in any direction)
      std::array<std::vector<double>,D> axis_points;
      for (int d = 0; d < D; d++) {
        int n = mesh.num_axis_points(d, Wonton::ALL);
        axis_points[d].resize(n);
        for (int i = 0; i < n; i++)
          axis_points[d][i] = mesh.get_axis_point(d, i-NG);
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

      // Check adjacencies
      verify_cell_adjacencies<D>(mesh_wrapper);
      verify_node_adjacencies<D>(mesh_wrapper);
    }
  }

} // namespace direct_product_mesh_wrapper_test 


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid1D) {
  const std::vector<double> axis_points1 = {0.0625, 0.125, 0.25, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<1>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid2D) {
  const std::vector<double> axis_points1 = {0.0, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<2>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, OneCell3D) {
  const std::vector<double> axis_points1 = {1.0, 1.5, 2.0};
  direct_product_mesh_wrapper_test::run_all_tests<3>(axis_points1);
}


// ============================================================================

TEST(Direct_Product_Mesh_Wrapper, SmallGrid10D) {
  const std::vector<double> axis_points1 = {0.0, 0.5, 1.0};
  direct_product_mesh_wrapper_test::run_all_tests<10>(axis_points1);
}


// ============================================================================

