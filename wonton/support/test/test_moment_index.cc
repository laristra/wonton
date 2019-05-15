/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include "wonton/support/moment_index.h"

// ============================================================================

template<int D>
void loop_over_exponents(
    int const order, int const d, std::array<int,D> & exponents, int & index) {
  int sum_exp =
    std::accumulate(exponents.begin(), exponents.begin()+d, 0);
  if (d == D-1) {
    // Base case
    // -- the last exponent is fixed by the need for sum(exponents) = order
    exponents[d] = order - sum_exp;
    // -- convert moment to index and confirm
    ASSERT_EQ(Wonton::moment_to_index<D>(order, exponents), index);
    // -- increment to next index
    index++;
  } else {
    // Recursive case
    // -- loop over all possible values for current exponent
    for (exponents[d] = order - sum_exp; exponents[d] >= 0; --exponents[d]) {
      // -- recurse to loop over next exponent
      loop_over_exponents<D>(order, d+1, exponents, index);
    }
  }
}

// ============================================================================
// Test core routine

template<int D>
void run_moment_index_test() {
  int order;
  std::array<int,D> exponents;
  // Confirm that certain hard-coded results are correct (index-to-moment)
  // -- index zero is always zeroth moment
  std::tie(order,exponents) = Wonton::index_to_moment<D>(0);
  ASSERT_EQ(order, 0);
  for (int d = 0; d < D; ++d) {
    ASSERT_EQ(exponents[d], 0);
  }
  // -- index 1..D are always first moments
  for (int d1 = 0; d1 < D; ++d1) {
    std::tie(order,exponents) = Wonton::index_to_moment<D>(d1+1);
    ASSERT_EQ(order, 1);
    for (int d2 = 0; d2 < D; ++d2) {
      if (d1 == d2) {
        ASSERT_EQ(exponents[d2], 1);
      } else {
        ASSERT_EQ(exponents[d2], 0);
      }
    }
  }
  // Confirm that certain hard-coded results are correct (moment-to-index)
  int index = 0;
  for (order = 0; order < 5; ++order) {
    loop_over_exponents<D>(order, 0, exponents, index);
  }
  // Confirm that results are reasonably self-consistent
  int const NMAX = 100;
  for (index = 0; index < NMAX; ++index) {
    // Convert index to moment specification
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

