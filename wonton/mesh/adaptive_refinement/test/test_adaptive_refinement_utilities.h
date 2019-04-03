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
#include "wonton/support/CellID.h"
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

// Getting tired of typing this
template<long D>
using BoxList = std::pair<std::vector<Wonton::CellID>,std::vector<Wonton::BoundingBox<D>>>;
//using BoxList = std::vector<std::pair<Wonton::CellID,Wonton::BoundingBox<D>>>;

// Getting tired of typing Wonton::
const int LO = Wonton::LO;
const int HI = Wonton::HI;

// ============================================================================
// Templates

template<long D>
double refinement_function(const Wonton::Point<D> r) {
  return(0.0);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<long D>
int num_cells() {
  return(0);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<long D>
BoxList<D> get_sample_points() {
  BoxList<D> boxes;
  return(std::move(boxes));
}

// ============================================================================
// 1D

template<>
double refinement_function<1>(const Wonton::Point<1> r) {
  return(2.0 + 3.0*r[0]);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
int num_cells<1>() {
  return(18);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
BoxList<1> get_sample_points<1>() {
  // Declare storage and counter
  std::vector<Wonton::CellID> id_list;
  std::vector<Wonton::BoundingBox<1>> box_list;
  int n;
  // First sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 0;
  box_list[n][0][LO] = 0.0;
  box_list[n][0][HI] = 1.0 / 8.0;
  // Second sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 6;
  box_list[n][0][LO] =  9.0 / 16.0;
  box_list[n][0][HI] = 10.0 / 16.0;
  // Third sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 14;
  box_list[n][0][LO] = 1.0 - 4.0 / 32.0;
  box_list[n][0][HI] = 1.0 - 3.0 / 32.0;
  // Return
  return(std::move(std::make_pair(id_list,box_list)));
}

// ============================================================================

template<>
double refinement_function<2>(const Wonton::Point<2> r) {
  return(1.2 - 1.1*r[0] + 2.4*r[1]);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
int num_cells<2>() {
  return(46);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
BoxList<2> get_sample_points<2>() {
  // Declare storage and counter
  std::vector<Wonton::CellID> id_list;
  std::vector<Wonton::BoundingBox<2>> box_list;
  int n;
  // First sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 4;
  box_list[n][0][LO] = 0.5;
  box_list[n][0][HI] = 1.0;
  box_list[n][1][LO] = 0.0;
  box_list[n][1][HI] = 0.5;
  // Second sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 20;
  box_list[n][0][LO] = 0.0;
  box_list[n][0][HI] = 1.0 / 16.0;
  box_list[n][1][LO] = 1.0 - 1.0 / 16.0;
  box_list[n][1][HI] = 1.0;
  // Third sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 35;
  box_list[n][0][LO] = 4.0 / 8.0;
  box_list[n][0][HI] = 5.0 / 8.0;
  box_list[n][1][LO] = 5.0 / 8.0;
  box_list[n][1][HI] = 6.0 / 8.0;
  // Return
  return(std::move(std::make_pair(id_list,box_list)));
}

// ============================================================================

template<>
double refinement_function<3>(const Wonton::Point<3> r) {
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

template<>
int num_cells<3>() {
  return(120);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
BoxList<3> get_sample_points<3>() {
  // Declare storage and counter
  std::vector<Wonton::CellID> id_list;
  std::vector<Wonton::BoundingBox<3>> box_list;
  int n;
  // First sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 0;
  box_list[n][0][LO] = 0.0;
  box_list[n][0][HI] = 0.25;
  box_list[n][1][LO] = 0.0;
  box_list[n][1][HI] = 0.25;
  box_list[n][2][LO] = 0.0;
  box_list[n][2][HI] = 0.25;
  // Second sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 92;
  box_list[n][0][LO] = 3.0 / 8.0;
  box_list[n][0][HI] = 4.0 / 8.0;
  box_list[n][1][LO] = 4.0 / 8.0;
  box_list[n][1][HI] = 5.0 / 8.0;
  box_list[n][2][LO] = 4.0 / 8.0;
  box_list[n][2][HI] = 5.0 / 8.0;
  // Third sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 119;
  box_list[n][0][LO] = 0.75;
  box_list[n][0][HI] = 1.0;
  box_list[n][1][LO] = 0.75;
  box_list[n][1][HI] = 1.0;
  box_list[n][2][LO] = 0.75;
  box_list[n][2][HI] = 1.0;
  // Return
  return(std::move(std::make_pair(id_list,box_list)));
}

// ============================================================================

template<>
double refinement_function<4>(const Wonton::Point<4> r) {
  return (2.0 - 3.0*r[0] + 1.2*r[1] - 0.5*r[2] + 0.3*r[3]);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
int num_cells<4>() {
  return(916);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

template<>
BoxList<4> get_sample_points<4>() {
  // Declare storage and counter
  std::vector<Wonton::CellID> id_list;
  std::vector<Wonton::BoundingBox<4>> box_list;
  int n;
  // First sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 0;
  box_list[n][0][LO] = 0.0;
  box_list[n][0][HI] = 0.25;
  box_list[n][1][LO] = 0.0;
  box_list[n][1][HI] = 0.25;
  box_list[n][2][LO] = 0.0;
  box_list[n][2][HI] = 0.25;
  box_list[n][3][LO] = 0.0;
  box_list[n][3][HI] = 0.25;
  // Second sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 92;
  box_list[n][0][LO] = 1.0 / 8.0;
  box_list[n][0][HI] = 2.0 / 8.0;
  box_list[n][1][LO] = 5.0 / 8.0;
  box_list[n][1][HI] = 6.0 / 8.0;
  box_list[n][2][LO] = 2.0 / 8.0;
  box_list[n][2][HI] = 3.0 / 8.0;
  box_list[n][3][LO] = 1.0 / 8.0;
  box_list[n][3][HI] = 2.0 / 8.0;
  // Third sample point
  // -- Expand lists
  n = id_list.size();
  id_list.resize(n+1);
  box_list.resize(n+1);
  // -- Save data
  id_list[n] = 119;
  box_list[n][0][LO] = 0.0;
  box_list[n][0][HI] = 1.0 / 8.0;
  box_list[n][1][LO] = 4.0 / 8.0;
  box_list[n][1][HI] = 5.0 / 8.0;
  box_list[n][2][LO] = 1.0 / 8.0;
  box_list[n][2][HI] = 2.0 / 8.0;
  box_list[n][3][LO] = 2.0 / 8.0;
  box_list[n][3][HI] = 3.0 / 8.0;
  // Return
  return(std::move(std::make_pair(id_list,box_list)));
}

// ============================================================================

}  // namespace Adaptive_Refinement_Utilities

#endif  // WONTON_ADAPTIVE_REFINEMENT_MESH_UTILITIES_H_
