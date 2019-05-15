/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_SUPPORT_COORDINATESYSTEM_H_
#define WONTON_SUPPORT_COORDINATESYSTEM_H_

#include <cmath>

#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

/*
  @file CoordinateSystem.h
  @brief Defines quantities needed for non-Cartesian coordinates

  Portage prefers to work with geometry factors of unity (actual volumes rather
  than wedges).  Because of this, a geometry factor is unnecessary.  However,
  to help emphasize to users of Portage that this is the assumption, the
  geometry factor is listed explicitly for every coordinate system.  This also
  makes it easier to change to different geometry factors in the future if
  desired.

  The modify_* routines all take an input calculated using the "standard",
  Cartesian-like expression, and then modify it to be appropriate for the
  selected coordinate system.  For Cartesian coordinates, this becomes a no-op
  because the expression is already calculated correctly.
 */


namespace Wonton {

  // Because C++ never bothered to define pi for some bizarre reason
  namespace CoordinateSystems {
    // Ensure we have pi to double precision
    double const minus_one = -1.0;
    double const pi = acos(minus_one);
    double const twopi = 2.0 * pi;
    double const fourpi = 2.0 * pi;
  }

// ============================================================================
/// Cartesian Coordinates
struct CartesianCoordinates {
 public:

  /// Geometry factor
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid for any positive dimensionality
    static_assert(D >= 1,
        "Cartesian coordinates must have positive dimensionality.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 0;

  /// Modify moments in-place to account for the coordinate system
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // No change from "standard", Cartesian-like calculation
    // --> Other than the geometry factor (which should be one, because any
    //     other value would be highly unusual for Cartesian coordinates).
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // No change from "standard", Cartesian-like calculation.
    // --> Other than the geometry factor (which should be one, because any
    //     other value would be highly unusual for Cartesian coordinates, but
    //     we verify this anyway).
    auto volume = vol0 * inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // No change from "standard", Cartesian-like calculation.
    // --> Other than the geometry factor (which should be one, because any
    //     other value would be highly unusual for Cartesian coordinates, but
    //     we verify this anyway).
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
  }

};  // Cartesian Coordinates


// ============================================================================
/// Cylindrical (Radial) Coordinates
/// Only valid in 1D.  The coordinate is the distance from the z axis.
struct CylindricalRadialCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 2 pi: a one-radian wedge of a cylinder.  Portage
  /// uses 1: a full cylinder.
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Cylindrical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 1;

  /// Modify moments in-place to account for the coordinate system
  ///   Cylindrical moments increment the radial coordinate exponent by one.
  /// That reduces the maximum order of the available moments by one.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // Find the maximum moment provided
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = index_to_moments<D>(moments.size()-1);
    // Given that we are reducing by a single order, how many moments will we
    // preserve?  In order words, how many moments are there from orders zero
    // to (max input order) - 1?
    int num = number_of_moments_through_order<D>(order-1);
    // Shift the moments down
    for (int n = 0; n < num; ++n) {
      // Get the moment specification
      std::tie(order,exponents) = index_to_moment<D>(n);
      // Increment
      order++;
      exponents[0]++;
      // Get the index
      int index = moment_to_index<D>(order, exponents);
      // Shift down
      moments[n] = moments[index];
    }
    // Resize array
    moments.resize(num);
    // Apply geometry factor
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = CoordinateSystems::twopi * mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    // -- Zeroth moment (volume)
    moments[0] *= CoordinateSystem::twopi * rhobar;
    // -- First moments (centroid * volume)
    moments[1] *= CoordinateSystem::twopi *
      (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Cylindrical (Radial) Coordinates


// ============================================================================
/// Cylindrical (Axisymmetric) Coordinates
/// Only valid in 2D.  The coordinates are the distance from the z axis and the
/// height.
struct CylindricalAxisymmetricCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 2 pi: a one-radian wedge of a cylinder.  Portage
  /// uses 1: a full cylinder.
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (axisymmetric) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 1;

  /// Modify moments in-place to account for the coordinate system
  ///   Cylindrical moments increment the radial coordinate exponent by one.
  /// That reduces the maximum order of the available moments by one.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // Find the maximum moment provided
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = index_to_moments<D>(moments.size()-1);
    // Given that we are reducing by a single order, how many moments will we
    // preserve?  In order words, how many moments are there from orders zero
    // to (max input order) - 1?
    int num = number_of_moments_through_order<D>(order-1);
    // Shift the moments down
    for (int n = 0; n < num; ++n) {
      // Get the moment specification
      std::tie(order,exponents) = index_to_moment<D>(n);
      // Increment
      order++;
      exponents[0]++;
      // Get the index
      int index = moment_to_index<D>(order, exponents);
      // Shift down
      moments[n] = moments[index];
    }
    // Resize array
    moments.resize(num);
    // Apply geometry factor
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = CoordinateSystems::twopi * mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    // -- Zeroth moment (volume)
    moments[0] *= CoordinateSystems::twopi * rhobar;
    // -- First moments (centroid * volume)
    moments[1] *= CoordinateSystems::twopi *
      (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[2] *= CoordinateSystems::twopi * rhobar;
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Cylindrical (Axisymmetric) Coordinates


// ============================================================================
/// Cylindrical (Polar) Coordinates
/// Only valid in 2D.  The coordinates are the distance from the z axis and the
/// azimuthal angle.
struct CylindricalPolarCoordinates {
 public:

  /// Geometry factor
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (polar) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 1;

  /// Modify moments in-place to account for the coordinate system
  ///   Cylindrical moments increment the radial coordinate exponent by one.
  /// That reduces the maximum order of the available moments by one.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // Find the maximum moment provided
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = index_to_moments<D>(moments.size()-1);
    // Given that we are reducing by a single order, how many moments will we
    // preserve?  In order words, how many moments are there from orders zero
    // to (max input order) - 1?
    int num = number_of_moments_through_order<D>(order-1);
    // Shift the moments down
    for (int n = 0; n < num; ++n) {
      // Get the moment specification
      std::tie(order,exponents) = index_to_moment<D>(n);
      // Increment
      order++;
      exponents[0]++;
      // Get the index
      int index = moment_to_index<D>(order, exponents);
      // Shift down
      moments[n] = moments[index];
    }
    // Resize array
    moments.resize(num);
    // Apply geometry factor
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    // -- Zeroth moment (volume)
    moments[0] *= rhobar;
    // -- First moments (centroid * volume)
    moments[1] *= (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[2] *= rhobar;
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Cylindrical Polar Coordinates


// ============================================================================
/// Cylindrical (3D) Coordinates
/// Only valid in 3D.  The coordinates are the distance from the z axis, the
/// azimuthal angle, and the height.
struct Cylindrical3DCoordinates {
 public:

  /// Geometry factor
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Cylindrical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 1;

  /// Modify moments in-place to account for the coordinate system
  ///   Cylindrical moments increment the radial coordinate exponent by one.
  /// That reduces the maximum order of the available moments by one.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // Find the maximum moment provided
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = index_to_moments<D>(moments.size()-1);
    // Given that we are reducing by a single order, how many moments will we
    // preserve?  In order words, how many moments are there from orders zero
    // to (max input order) - 1?
    int num = number_of_moments_through_order<D>(order-1);
    // Shift the moments down
    for (int n = 0; n < num; ++n) {
      // Get the moment specification
      std::tie(order,exponents) = index_to_moment<D>(n);
      // Increment
      order++;
      exponents[0]++;
      // Get the index
      int index = moment_to_index<D>(order, exponents);
      // Shift down
      moments[n] = moments[index];
    }
    // Resize array
    moments.resize(num);
    // Apply geometry factor
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    // -- Zeroth moment (volume)
    moments[0] *= rhobar;
    // -- First moments (centroid * volume)
    moments[1] *= (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[2] *= rhobar;
    moments[3] *= rhobar;
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Cylindrical (3D) Coordinates


// ============================================================================
/// Spherical (Radial) Coordinates
/// Only valid in 1D.  The coordinate is the distance from the origin.
struct SphericalRadialCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 4 pi: a one-steradian wedge of a sphere.  Portage
  /// uses 1: a full sphere.
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Spherical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 2;

  /// Modify moments in-place to account for the coordinate system
  ///   Spherical moments increment the radial coordinate exponent by two.
  /// That reduces the maximum order of the available moments by two.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    // Find the maximum moment provided
    int order;
    std::array<int,D> exponents;
    std::tie(order,exponents) = index_to_moments<D>(moments.size()-1);
    // Given that we are reducing by a single order, how many moments will we
    // preserve?  In order words, how many moments are there from orders zero
    // to (max input order) - 1?
    int num = number_of_moments_through_order<D>(order-1);
    // Shift the moments down
    for (int n = 0; n < num; ++n) {
      // Get the moment specification
      std::tie(order,exponents) = index_to_moment<D>(n);
      // Increment
      order += 2;
      exponents[0] += 2;
      // Get the index
      int index = moment_to_index<D>(order, exponents);
      // Shift down
      moments[n] = moments[index];
    }
    // Resize array
    moments.resize(num);
    // Apply geometry factor
    for (int n = 0; n < moments.size(); ++n) {
      moments[n] *= inv_geo_fac;
    }
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (plo[0] + phi[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    auto volume = CoordinateSystems::fourpi *
      (rbar*rbar + dr_2*dr_2/3.0) * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (phi[0] + plo[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    // -- Zeroth moment (volume)
    moments[0] *= CoordinateSystems::fourpi * (rbar*rbar + dr_2*dr_2/3.0);
    // -- First moments (centroid * volume)
    moments[1] *= CoordinateSystems::fourpi * (rbar*rbar + dr_2*dr_2);
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Spherical (Radial) Coordinates


// ============================================================================
/// Spherical (3D) Coordinates
/// Only valid in 3D.  The coordinates are the distance from the origin, the
/// angle of inclination, and the azimuthal angle.
struct Spherical3DCoordinates {
 public:

  /// Geometry factor
  static constexpr double geometry_factor = 1;

  /// Inverse geometry factor
  static constexpr double inv_geo_fac = 1.0 / geometry_factor;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Spherical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static void modify_gradient(Wonton::Vector<D> & gradient,
      Wonton::Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
    gradient[2] /= (reference_point[0] * sin(reference_point[1]));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static void modify_line_element(Wonton::Vector<D> & line_element,
      Wonton::Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
    line_element[2] *= (reference_point[0] * sin(reference_point[1]));
  }

  /// How many orders of moments do you lose?
  static constexpr int moment_shift = 0;

  /// Modify moments in-place to account for the coordinate system
  ///   You cannot use the moment-shift method with 3D spherical coordinates,
  /// because the volume element includes a sin(theta) term.
  template<long D>
  static void modify_moments(std::vector<double> & moments) {
    static_assert(false,
        "Moment-shifting cannot be used with 3D spherical coordinates.")
  }

  /*
  /// Modify volume to account for the coordinate system
  template<long D>
  static double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (plo[0] + phi[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    thetabar = 0.5 * (phi[1] + plo[1]);
    dtheta_2 = 0.5 * (phi[1] - plo[1]);
    sin_tb  = sin(thetabar);
    sinc_dt = sin(dtheta_2) / dtheta_2;
    auto volume *= (rbar*rbar + dr_2*dr_2/3.0) * (sin_tb * sinc_dt) * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return(volume);
  }
  */

  /// Modify moments to account for the coordinate system
  ///   This is an optimization for cells that are axis-aligned boxes, in order
  /// to avoid computing higher-order moments, then throwing away most or all
  /// of them.  This routine is hard-coded to assume the zeroth and first
  /// moments, because it is not obvious that higher-order moments are needed
  /// and (for this optimization) it is necessary to explicitly derive the
  /// appropriate expressions.
  template<long D>
  static void modify_moments(std::vector<double> const & moments,
      Wonton::Point<D> const & plo, Wonton::Point<D> const phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (phi[0] + plo[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    rbar_sq = rbar * rbar;
    dr_2_sq = dr_2 * dr_2;
    rr1 = rbar_sq + dr_2_sq;
    rr2 = rbar_sq + dr_2_sq/3.0;
    thetabar = 0.5 * (phi[1] + plo[1]);
    dtheta_2 = 0.5 * (phi[1] - plo[1]);
    sin_tb  = sin(thetabar);
    sinc_dt = sin(dtheta_2) / dtheta_2;
    cosc_tb = cos(thetabar) / thetabar;
    cos_dt  = cos(dtheta_2);
    ss1 = sin_tb * sinc_dt;
    // -- Zeroth moment (volume)
    moments[0] *= (rbar*rbar + dr_2*dr_2/3.0) * (sin_tb * sinc_dt);
    // -- First moments (centroid * volume)
    moments[1] *= rr1 * ss1;
    moments[2] *= rr2 * (ss1 + cosc_tb * (sinc_dt - cos_dt));
    moments[3] *= rr2 * ss1;
    // Apply geometry factor
    for (int d = 0; d < 1+D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return(std::move(moments));
  }

};  // Spherical (3D) Coordinates


// ============================================================================

// Default coordinate system for consistency across the code
using DefaultCoordSys = CartesianCoordinates;

// ============================================================================

}  // namespace Wonton

#endif  // WONTON_SUPPORT_COORDINATESYSTEM_H_
