/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_ADAPTIVE_REFINEMENT_MESH_UTILITIES_H_
#define WONTON_ADAPTIVE_REFINEMENT_MESH_UTILITIES_H_

#include <iostream>
#include <cmath>
#include <vector>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

#include "gtest/gtest.h"

// ============================================================================
/*!
  @file test_adaptive_refinement_utilities.cc
  @brief Support function for testing the Adaptive_Refinement_Mesh and its
         wrapper.

  NOTE 1: The tests are calibrated against these refinement functions.  If you
  change these functions, you will have to recalibrate the tests.
*/

namespace Adaptive_Refinement_Utilities {

// ============================================================================
// 1D

template<long D>
double refinement_function(const Wonton::Point<1> r) {
  return(2.0 + 3.0*r[0]);
}

// ============================================================================

template<long D>
double refinement_function(const Wonton::Point<2> r) {
  return(1.2 - 1.1*r[0] + 2.4*r[1]);
}

// ============================================================================

template<long D>
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

// ============================================================================

template<long D>
double refinement_function(const Wonton::Point<4> r) {
  return (2.0 - 3.0*r[0] + 1.2*r[1] - 0.5*r[2] + 0.3*r[3]);
}

// ============================================================================

}  // namespace Adaptive_Refinement_Utilities

#endif  // WONTON_ADAPTIVE_REFINEMENT_MESH_UTILITIES_H_
