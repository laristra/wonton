/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <numeric>

#include "gtest/gtest.h"

#include "wonton/support/moment_index.h"

// ============================================================================
// Data for hard-coded examples

template<int D>
constexpr int get_index() {
  return 17;
}

// ----------------------------------------------------------------------------

template<int D>
constexpr int get_order();

template<>
constexpr int get_order<1>() {
  return 17;
}

template<>
constexpr int get_order<2>() {
  return 5;
}

template<>
constexpr int get_order<3>() {
  return 3;
}

// ----------------------------------------------------------------------------

template<int D>
constexpr std::array<int,D> get_exponents();

template<>
constexpr std::array<int,1> get_exponents<1>() {
  std::array<int,1> exponents{17};
  return exponents;
}

template<>
constexpr std::array<int,2> get_exponents<2>() {
  std::array<int,2> exponents{3,2};
  return exponents;
}

template<>
constexpr std::array<int,3> get_exponents<3>() {
  std::array<int,3> exponents{0,2,1};
  return exponents;
}

// ============================================================================
// Test core routine

template<int D>
void run_moment_index_test() {
  int const NMAX = 100;
  for (int index = 0; index < NMAX; ++index) {
    // Convert index to moment specification
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = Wonton::index_to_moment<D>(index);
    // Verify moment specification is internally consistent
    ASSERT_GE(order, 0);
    for (int d = 0; d < D; ++d) {
      ASSERT_GE(exponents[d], 0);
    }
    ASSERT_EQ(order, std::accumulate(exponents.begin(), exponents.end(), 0));
    // Convert moment specification back to index
    auto index2 = Wonton::moment_to_index<D>(order, exponents);
    // Ensure the round trip is self-consistent
    ASSERT_EQ(index, index2);
  }

  // Test sample point(s)
  int index_in = get_index<D>();
  int order_in = get_order<D>();
  std::array<int,D> exponents_in = get_exponents<D>();
  int order_out;
  std::array<int,D> exponents_out;
  std::tie(order_out,exponents_out) = Wonton::index_to_moment<D>(index_in);
  int index_out = Wonton::moment_to_index<D>(order_in, exponents_in);
  ASSERT_EQ(index_in, index_out);
  ASSERT_EQ(order_out, order_out);
  for (int d = 0; d < D; ++d) {
    ASSERT_EQ(exponents_in[d], exponents_out[d]);
  }
  int order0;
  std::array<int,D> exponents0;
  std::tie(order0,exponents0) = Wonton::index_to_moment<D>(0);
  ASSERT_EQ(order0,0);
  for (int d = 0; d < D; ++d) {
    ASSERT_EQ(exponents0[d],0);
  }
}

// ============================================================================

// 1D
TEST(Moment_Index_test, Test1D) {
  run_moment_index_test<1>();
}

// 2D
TEST(Moment_Index_test, Test2D) {
  run_moment_index_test<2>();
}

// 3D
TEST(Moment_Index_test, Test3D) {
  run_moment_index_test<3>();
}

