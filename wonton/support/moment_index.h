/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/
#ifndef WONTON_SUPPORT_MOMENT_INDEX_
#define WONTON_SUPPORT_MOMENT_INDEX_

#include <array>
#include <cassert>
#include <cmath>
#include <utility>

//! This file contains convenience functions to help extract moments from the
//! moments list by converting between a moment specification and a moments
//! list index.
//!
//! A moment specification consists of two parts: an order and the exponents of
//! the basis axes.  Consider some 2D Cartesian examples, so that the basis
//! axes are x and y.  The sum of the exponents is required to be equal to the
//! order.  So if, for example, you are asking for the moment
//! integral{x y^2 dx dy} then the moment specification would be order = 3 and
//! exponents = {1, 2}.  If you wanted to enumerate all 2nd-order moments, then
//! the options would be:
//!    order = 2, exponents = {0, 2}
//!    order = 2, exponents = {1, 1}
//!    order = 2, exponents = {2, 0}
//!
//! A more complex example could be 3D cylindrical coordinates.  Because this
//! is in cylindrical coordinates, the basis axes are rho, phi, and z.  The
//! second-order moments would be:
//!    order = 2, exponents = {0, 0, 2}
//!    order = 2, exponents = {0, 2, 0}
//!    order = 2, exponents = {2, 0, 0}
//!    order = 2, exponents = {0, 1, 1}
//!    order = 2, exponents = {1, 0, 1}
//!    order = 2, exponents = {1, 1, 0}
//! However, because the volume element dV is equal to drho rho dphi dz, then
//! the moment with exponents {0,1,1} is be computed as
//! integral{phi z dV} = integral{rho phi z drho dphi dz} rather than as
//! integral{phi z drho dphi dz}.

namespace Wonton {

  // ==========================================================================

  //! Convert from moment specification to index
  //! @param[in]  order  Order of the moment
  //! @param[in]  exponents  Powers for the coordinate axes
  template<int D>
  constexpr int moment_to_index(
      int const order, std::array<int,D> const exponents);

  template<>
  constexpr int moment_to_index<3>(
      int const order, std::array<int,3> const exponents) {
    int constexpr D = 3;
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    int index = 0;
    index += (order+2) * (order+1) * order / 6;
    int order1 = order - exponents[0];
    index += (order1+1) * order1 / 2;
    index += order1 - exponents[1];
    return(index);
  }

  template<>
  constexpr int moment_to_index<2>(
      int const order, std::array<int,2> const exponents) {
    int constexpr D = 2;
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    int index = 0;
    index += (order+1) * order / 2;
    index += order - exponents[0];
    return(index);
  }

  template<>
  constexpr int moment_to_index<1>(
      int const order, std::array<int,1> const exponents) {
    int constexpr D = 1;
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    int index = 0;
    index += order;
    return(index);
  }

  // ==========================================================================

  //! Count the number of moments up to specified order
  //! @param[in] order  The order up to which to count the moments (inclusive)
  template<int D>
  constexpr int count_moments(int const order);

  template<>
  constexpr int count_moments<3>(int const order) {
    return (order+3) * (order+2) * (order+1) / 6;
  }

  template<>
  constexpr int count_moments<2>(int const order) {
    return (order+2) * (order+1) / 2;
  }

  template<>
  constexpr int count_moments<1>(int const order) {
    return order+1;
  }

  // ==========================================================================
  // TODO: We can't get a non-const reference to a std::array element in a
  //       constexpr function until C++17.  Thus these routines cannot be
  //       constexpr until Wonton starts using C++17.  Once these become
  //       constexpr the inline specifier is redundant.

  //! Convert from index to moment specification
  //! @param[in]  index  Index in moment list
  template<int D>
  inline std::pair<int,std::array<int,D>> index_to_moment(int const index);

  template<>
  inline std::pair<int,std::array<int,3>> index_to_moment<3>(
     int const index) {
    int constexpr D = 3;
    int n = index;
    // Declare the output variables
    int order = 0;
    std::array<int,D> exponents = {};
    // Find order
    //    While this could be done analytically, the expression is ugly and
    // involves cube roots and complex numbers.  With only a small number of
    // orders, this should be plenty fast enough (and possibly faster than the
    // ugly expression).
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
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    // Return
    return std::make_pair(order,exponents);
  }

  template<>
  inline std::pair<int,std::array<int,2>> index_to_moment<2>(
      int const index) {
    int constexpr D = 2;
    int n = index;
    // Declare the output variables
    int order = 0;
    std::array<int,D> exponents = {};
    // Find order
    order = (int) floor((sqrt((double)8*n+1) - 1) / 2);
    // Skip the lower-order moments
    n -= (order+1) * order / 2;
    // Find the first exponent
    exponents[0] = order - n;
    // Find the second exponent
    exponents[1] = n;
    // Verify results
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    // Return
    return std::make_pair(order,exponents);
  }

  template<>
  inline std::pair<int,std::array<int,1>> index_to_moment<1>(
      int const index) {
    int constexpr D = 1;
    // Declare the output variables
    int order = 0;
    std::array<int,D> exponents = {};
    // Find order
    order = index;
    // Find the first exponent
    exponents[0] = index;
    // Verify results
    int exponent_sum = 0;
    for (int d = 0; d < D; ++d) {
      assert(exponents[d] >= 0);
      exponent_sum += exponents[d];
    }
    assert(exponent_sum == order);
    // Return
    return std::make_pair(order,exponents);
  }

  // ==========================================================================
  // TODO: Experimental recursive implementations below.

  /*
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

  // Specialization for base case
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

  // Specialization for base case
  template<>
  int get_moment_order_from_index<2>(int const index) {
    double radicand = 8*index + 1;
    double order = floor((sqrt(radicand) - 1) / 2);
    return((int) order);
  }

  // Specialization for base case
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

  // Specialization for base case
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

  // Specialization for base case
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
    return std::make_pair(order,exponents);
  }
  */

}  // namespace Wonton

#endif  // WONTON_SUPPORT_MOMENT_INDEX_
