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

namespace Wonton {

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
    return(std::move(std::make_pair(order,exponents)));
  }

}  // namespace Wonton

#endif  // PORTAGE_SUPPORT_MOMENT_INDEX_
