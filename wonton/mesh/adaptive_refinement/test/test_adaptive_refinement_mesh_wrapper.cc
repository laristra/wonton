/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh_wrapper.h"
#include "wonton/mesh/adaptive_refinement/test/test_adaptive_refinement_utilities.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

// ============================================================================
/*!
  @file test_adaptive_refinement_mesh_wrapper.cc
  @brief Tests for the adaptive_refinement mesh wrapper definition in
         adaptive_refinement_mesh_wrapper.h
*/

namespace adaptive_refinement_mesh_test {

template <int D>
void run_basic_tests() {
  // Create a mesh
  // -- The tests are calibrated to a level-zero mesh ranging from 0 to 1 along
  //    every axis.
  Wonton::Point<D> lo, hi;
  for (int d = 0; d < D; ++d) {
    lo[d] = 0.0;
    hi[d] = 1.0;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(
      &Adaptive_Refinement_Utilities::refinement_function<D>, lo, hi);

  // Create a wrapper
  Wonton::Adaptive_Refinement_Mesh_Wrapper<D> wrapper(mesh);

  // Dimensionality
  ASSERT_EQ(wrapper.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(wrapper.total_num_cells(),
      Adaptive_Refinement_Utilities::num_cells<D>());

  // Cell coordinates
  // -- These are known from testing
  Adaptive_Refinement_Utilities::BoxList<D> boxes =
    Adaptive_Refinement_Utilities::get_sample_points<D>();
  auto id_list = std::get<std::vector<Wonton::CellID>>(boxes);
  auto box_list = std::get<std::vector<Wonton::BoundingBox<D>>>(boxes);
  ASSERT_TRUE(id_list.size() > 0); // Ensure that we got the explicit
                                   // specialization, rather than an automatic
                                   // specialization of the template.
  ASSERT_EQ(id_list.size(), box_list.size());
  for (int n = 0; n < id_list.size(); ++n) {
    Wonton::Point<D> plo, phi;
    wrapper.cell_get_bounds(id_list[n], &plo, &phi);
    auto expected = box_list[n];
    for (int d = 0; d < D; ++d) {
      ASSERT_EQ(plo[d], expected[d][Wonton::LO]);
      ASSERT_EQ(phi[d], expected[d][Wonton::HI]);
    }
  }
}   // void run_basic_tests

}


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test1D) {
  adaptive_refinement_mesh_test::run_basic_tests<1>();
}

// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test2D) {
  adaptive_refinement_mesh_test::run_basic_tests<2>();
}

// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test3D) {
  adaptive_refinement_mesh_test::run_basic_tests<3>();
}

// ============================================================================
// Not sure why you'd want a 4D AMR-like mesh, but this proves that it works
// (verifying that I didn't hard-code any assumptions about the dimensionality
// being <= 3, which I was worried about).

TEST(Adaptive_Refinement_Mesh, Test4D) {
  adaptive_refinement_mesh_test::run_basic_tests<4>();
}

