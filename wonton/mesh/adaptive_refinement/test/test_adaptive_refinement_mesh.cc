/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"

#include <iostream>
#include <cmath>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

// ============================================================================
/*!
  @file test_adaptive_refinement_mesh.cc
  @brief Tests for the adaptive_refinement mesh definition in
         adaptive_refinement_mesh.h
*/

// ============================================================================
// NOTE: The tests below are calibrated against these refinement functions.  If
//       you change these functions, you will have to recalibrate the tests.
// NOTE: These all have the same name because the templating of the grid on
//       different dimensionalities provides enough information for the
//       compiler to select the correct function.

double refinement_function(const Wonton::Point<1> r) {
  return(2.0 + 3.0*r[0]);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double refinement_function(const Wonton::Point<2> r) {
  return(1.2 - 1.1*r[0] + 2.4*r[1]);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double refinement_function(const Wonton::Point<3> r) {
  const int DIM = 3;
  double radius = 0;
  for (int d = 0; d < DIM; ++d) {
    double x = r[d] - 0.5;
    radius += x * x;
  }
  radius = sqrt(radius);
  return(2.5 - 6.0*radius*radius);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

double refinement_function(const Wonton::Point<4> r) {
  return (2.0 - 3.0*r[0] + 1.2*r[1] - 0.5*r[2] + 0.3*r[3]);
}

// ============================================================================
// Not sure why you'd want a 4D AMR-like mesh, but this proves that it works
// (verifying that I didn't hard-code any assumptions about the dimensionality
// being <= 3, which I was worried about).

TEST(Adaptive_Refinement_Mesh, Test4D) {

  const int D = 4;

  // Create a single cell mesh
  Wonton::Point<D> lo, hi;
  for (int d = 0; d < D; ++d) {
    lo[d] = 0.0;
    hi[d] = 1.0;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 916);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(0);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 0.25);
  ASSERT_EQ(bounds[1][0], 0.0);
  ASSERT_EQ(bounds[1][1], 0.25);
  ASSERT_EQ(bounds[2][0], 0.0);
  ASSERT_EQ(bounds[2][1], 0.25);
  ASSERT_EQ(bounds[3][0], 0.0);
  ASSERT_EQ(bounds[3][1], 0.25);
  bounds = mesh.cell_get_bounds(92);
  ASSERT_EQ(bounds[0][0], 1.0/8.0);
  ASSERT_EQ(bounds[0][1], 1.0/4.0);
  ASSERT_EQ(bounds[1][0], 5.0/8.0);
  ASSERT_EQ(bounds[1][1], 3.0/4.0);
  ASSERT_EQ(bounds[2][0], 1.0/4.0);
  ASSERT_EQ(bounds[2][1], 3.0/8.0);
  ASSERT_EQ(bounds[3][0], 1.0/8.0);
  ASSERT_EQ(bounds[3][1], 1.0/4.0);
  bounds = mesh.cell_get_bounds(119);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 1.0/8.0);
  ASSERT_EQ(bounds[1][0], 1.0/2.0);
  ASSERT_EQ(bounds[1][1], 5.0/8.0);
  ASSERT_EQ(bounds[2][0], 1.0/8.0);
  ASSERT_EQ(bounds[2][1], 1.0/4.0);
  ASSERT_EQ(bounds[3][0], 1.0/4.0);
  ASSERT_EQ(bounds[3][1], 3.0/8.0);

}


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test3D) {

  const int D = 3;

  // Create a single cell mesh
  Wonton::Point<D> lo, hi;
  for (int d = 0; d < D; ++d) {
    lo[d] = 0.0;
    hi[d] = 1.0;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 120);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(0);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 0.25);
  ASSERT_EQ(bounds[1][0], 0.0);
  ASSERT_EQ(bounds[1][1], 0.25);
  ASSERT_EQ(bounds[2][0], 0.0);
  ASSERT_EQ(bounds[2][1], 0.25);
  bounds = mesh.cell_get_bounds(92);
  ASSERT_EQ(bounds[0][0], 3.0/8.0);
  ASSERT_EQ(bounds[0][1], 1.0/2.0);
  ASSERT_EQ(bounds[1][0], 1.0/2.0);
  ASSERT_EQ(bounds[1][1], 5.0/8.0);
  ASSERT_EQ(bounds[2][0], 1.0/2.0);
  ASSERT_EQ(bounds[2][1], 5.0/8.0);
  bounds = mesh.cell_get_bounds(119);
  ASSERT_EQ(bounds[0][0], 0.75);
  ASSERT_EQ(bounds[0][1], 1.0);
  ASSERT_EQ(bounds[1][0], 0.75);
  ASSERT_EQ(bounds[1][1], 1.0);
  ASSERT_EQ(bounds[2][0], 0.75);
  ASSERT_EQ(bounds[2][1], 1.0);

}


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test2D) {

  const int D = 2;

  // Create a mesh
  Wonton::Point<D> lo, hi;
  for (int d = 0; d < D; ++d) {
    lo[d] = 0.0;
    hi[d] = 1.0;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 46);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(4);
  ASSERT_EQ(bounds[0][0], 0.5);
  ASSERT_EQ(bounds[0][1], 1.0);
  ASSERT_EQ(bounds[1][0], 0.0);
  ASSERT_EQ(bounds[1][1], 0.5);
  bounds = mesh.cell_get_bounds(20);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 1.0/16.0);
  ASSERT_EQ(bounds[1][0], 1.0-1.0/16.0);
  ASSERT_EQ(bounds[1][1], 1.0);
  bounds = mesh.cell_get_bounds(35);
  ASSERT_EQ(bounds[0][0], 0.5);
  ASSERT_EQ(bounds[0][1], 5.0/8.0);
  ASSERT_EQ(bounds[1][0], 5.0/8.0);
  ASSERT_EQ(bounds[1][1], 3.0/4.0);

}


// ============================================================================

TEST(Adaptive_Refinement_Mesh, Test1D) {

  const int D = 1;

  // Create a mesh
  Wonton::Point<D> lo, hi;
  for (int d = 0; d < D; ++d) {
    lo[d] = 0.0;
    hi[d] = 1.0;
  }
  Wonton::Adaptive_Refinement_Mesh<D> mesh(&refinement_function, lo, hi);

  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), D);

  // Cell counts
  // -- This is known from testing
  ASSERT_EQ(mesh.num_cells(), 18);

  // Cell coordinates
  // -- These are known from testing
  auto bounds = mesh.cell_get_bounds(0);
  ASSERT_EQ(bounds[0][0], 0.0);
  ASSERT_EQ(bounds[0][1], 0.125);
  bounds = mesh.cell_get_bounds(6);
  ASSERT_EQ(bounds[0][0], 9.0/16.0);
  ASSERT_EQ(bounds[0][1], 10.0/16.0);
  bounds = mesh.cell_get_bounds(14);
  ASSERT_EQ(bounds[0][0], 1.0 - 4.0/32.0);
  ASSERT_EQ(bounds[0][1], 1.0 - 3.0/32.0);

}


