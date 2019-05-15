/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef PORTAGE_SUPPORT_MOMENT_INDEX_
#define PORTAGE_SUPPORT_MOMENT_INDEX_

#include <array>
#include <cassert>
#include <cmath>
#include <utility>

//! This file contains convenience functions to help extract moments from the
//! moments list by converting between a moment specification and a moments
//! list index.

// -- index in moment list for moment of order M with D exponents exp: index
//    index<D> = sum[d1=1..D](N<d1>(M-1-sum[d2=1..d-1](exp[d2])))
// It looks as through the index function could be defined recursively:
//    int moment_to_index<D>(order, exponents) {
//      auto term1 = number_of_moments_through_order<D>(order - 1);
//      auto term2 = moment_to_index<D>(order-exp[0], exp[1:])
//      return(term1 + term2);
//    }
// This could might be able to be simplified using a variadic template
// recursive function, but that might make the original call difficult.  Then
// you need the base case:
//    int moment_to_index<0>(order, exponents) {
//      return(0);
//    }
// If this is true, then it would be possible to just write a single recursive
// template function for number_of_moments_through_order<D>,
// moment_to_index<D>, and perhaps even index_to_moment<D> (at least as a
// search if not as an invertible function).  These inversion functions could
// also be stand-alone as inverses of number_of_moments_through_order and be
// used to figure out the maximum order of a list without having to construct
// the entire moment specification, but stopping after figuring out the order
// only.

namespace Portage {

  // ==========================================================================

  //! Number of moments for orders zero through specified order
  //! @param[in]  o_max  Maximum order
  // TODO: This implementation is not yet confirmed to be correct for D > 3.  I
  //       think it's right, but at the moment it's more of an educated guess.
  template<int D>
  int number_of_moments_through_order(int const o_max) {
    static_assert(D > 0,
        "\"number_of_moments_through_order\" only works for D > 0");
    // Ordering could be important because of integer truncation
    auto term_d1 = number_of_moments_through_order<D-1>(o_max);
    return(term_d1 * (o_max + D) / D);
  }

  template<>
  int number_of_moments_through_order<1>(int const o_max) {
    return(o_max+1);
  }

  // ==========================================================================

  //! Get moment order from index
  //! @param[in]  index  Index in moment list
  // TODO: This implementation is not yet confirmed to be correct for D > 3.  I
  //       think it's right, but at the moment it's more of an educated guess.
  template<int D>
  int get_moment_order_from_index(int const index) {
    static_assert(D > 0,
        "\"number_of_moments_through_order\" only works for D > 0");
    int order = 0;
    while(true) {
      if (index < number_of_moments_through_order<D>(order))
        break;
      order++;
    }
    return(order);
  }

  template<>
  int get_moment_order_from_index<2>(int const index) {
    double radicand = 8*index + 1;
    double order = floor((sqrt(radicand) - 1) / 2);
    return((int) order);
  }

  template<>
  int get_moment_order_from_index<1>(int const index) {
    return(index);
  }

  // ==========================================================================

  // Utility routine to do recursive calculation
  // TODO: This implementation is not yet confirmed to be correct for D > 3.  I
  //       think it's right, but at the moment it's more of an educated guess.
  template<int D>
  int moment_to_index_recursive(
      int const order, int const * const exponents) {
    auto term1 = number_of_moments_through_order<D>(order-1);
    auto term2 = moment_to_index_recursive<D-1>(
        order-(*exponents), exponents+1);
    return(term1 + term2);
  }

  template<>
  int moment_to_index_recursive<1>(
      int const order, int const * const exponents) {
    return(number_of_moments_through_order<1>(order-1));
  }

  //! Convert from moment specification to index
  //! @param[in]  order  Order of the moment
  //! @param[in]  exponents  Powers for the coordinate axes
  template<int D>
  int moment_to_index(int const order, std::array<int,D> const exponents) {
    static_assert(D > 0,
        "\"moment_to_index_recursive\" only works for D > 0");
    return(moment_to_index_recursive<D>(order, exponents.data()));
  }
/*
  template<>
  int moment_to_index<3>(int const order, std::array<int,3> const exponents) {
    int constexpr D = 3;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += (order+2) * (order+1) * order / 6;
    int order1 = order - exponents[0];
    index += (order1+1) * order1 / 2;
    index += order1 - exponents[1];
    return(index);
  }

  template<>
  int moment_to_index<2>(int const order, std::array<int,2> const exponents) {
    int constexpr D = 2;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += (order+1) * order / 2;
    index += order - exponents[0];
    return(index);
  }

  template<>
  int moment_to_index<1>(int const order, std::array<int,1> const exponents) {
    int constexpr D = 1;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    int index = 0;
    index += order;
    return(index);
  }
*/
  // ==========================================================================

  // Utility routine to do recursive calculation
  // TODO: This implementation is not yet confirmed to be correct for D > 3.  I
  //       think it's right, but at the moment it's more of an educated guess.
  template<int D>
  void index_to_moment_recursive(
      int const index, int & order, int * exponents) {
    order = get_moment_order_from_index<D>(index);
    int reduced_index = index - number_of_moments_through_order<D>(order-1);
    int temp; // order - exponents[0];
    index_to_moment_recursive<D-1>(reduced_index, temp, exponents+1);
    (*exponents) = order - temp;
  }
  /*
  Example: M = 3, i = 1, j = 0, k = 2  -->  n = 15
  <3> : index = 15
        order = 3
        N<3>(3-1) = 10
        reduced_index = 5
        <2> : index = 5
              order = 2
              N<2>(2-1) = 3
              reduced_index = 2
              <1> : index = 2
                    order = 2 = M - i - j
                    exponents = [2] = [k]
              temp = 2, exponents = [#,2]
              exponents = [0,2]
        temp = 2, exponents = [#,0,2]
        exponents = [1,0,2]
  index = 15, order = 3, exponents = [1,0,2] = [i,j,k]
  */

  template<>
  void index_to_moment_recursive<1>(
      int const index, int & order, int * exponents) {
    order = index;
    (*exponents) = index;
  }

  //! Convert from index to moment specification
  //! @param[in]  index  Index in moment list
  template<int D>
  std::pair<int,std::array<int,D>> index_to_moment(int const index) {
    static_assert(D > 0,
        "\"moment_to_index_recursive\" only works for D > 0");
    int order;
    std::array<int,D> exponents;
    index_to_moment_recursive<D>(index, order, exponents.data());
    return(std::move(std::make_pair(order,exponents)));
  }

/*
  //! Convert from index to moment specification
  //! @param[in]  index  Index in moment list
  template<int D>
  std::pair<int,std::array<int,D>> index_to_moment(int const index) {
    static_assert(false, "\"index_to_moment\" only works in 1, 2, or 3D");
  }

  template<>
  std::pair<int,std::array<int,3>> index_to_moment<3>(int const index) {
    int constexpr D = 3;
    int n = index;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    //    While this could be done analytically, the expression is ugly and
    // involves cube roots and complex numbers.  With only a small number of
    // orders, this should be plenty fast enough (and possibly faster than the
    // ugly expression).
    order = 0;
    while (true) {
      if (n < (order+2)*(order+1)*order/6)
        break;
      order++;
    }
    order--;
    // Skip the lower-order moments
    n -= (order+2)*(order+1)*order/6;
    // Find the first exponent
    exponents[0] = order - (int) floor((sqrt((double)8*n+1) - 1) / 2);
    auto order1 = order - exponents[0];
    // Skip moments with wrong first exponent
    n -= (order1+1) * order1 / 2;
    // Find the second exponent
    exponents[1] = order1 - n;
    // Find the third exponent
    exponents[2] = n;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }

  template<>
  std::pair<int,std::array<int,2>> index_to_moment<2>(int const index) {
    int constexpr D = 2;
    int n = index;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    order = (int) floor((sqrt((double)8*n+1) - 1) / 2);
    // Skip the lower-order moments
    n -= (order+1) * order / 2;
    // Find the first exponent
    exponents[0] = order - n;
    // Find the second exponent
    exponents[1] = n;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }

  template<>
  std::pair<int,std::array<int,1>> index_to_moment<1>(int const index) {
    int constexpr D = 1;
    // Declare the output variables
    int order;
    std::array<int,D> exponents;
    // Find order
    order = index;
    // Find the first exponent
    exponents[0] = index;
    // Verify results
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
    }
    auto exponent_sum = std::accumulate(exponents.begin(), exponents.end(), 0);
    assert(exponent_sum == order);
    // Return
    return(std::move(std::make_pair(order,exponents)));
  }
*/

}  // namespace Portage

#endif  // PORTAGE_SUPPORT_MOMENT_INDEX_
