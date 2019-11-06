/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_SUPPORT_COORDINATESYSTEM_H_
#define WONTON_SUPPORT_COORDINATESYSTEM_H_

#include <cmath>

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
    // atan, acos, and similar are not (by the standard) constexpr, so we can't
    // use them to define constexpr values for pi
    constexpr double pi = 3.141592653589793238462643383279502884L;
    constexpr double twopi = 2.0 * pi;
    constexpr double fourpi = 2.0 * pi;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid for any positive dimensionality
    static_assert(D >= 1,
        "Cartesian coordinates must have positive dimensionality.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // No change from "standard", Cartesian-like calculation.
    // --> Other than the geometry factor (which should be one, because any
    //     other value would be highly unusual for Cartesian coordinates, but
    //     we verify this anyway).
    auto volume = vol0 * inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // No change from "standard", Cartesian-like calculation.
    // --> Other than the geometry factor (which should be one, because any
    //     other value would be highly unusual for Cartesian coordinates, but
    //     we verify this anyway).
    auto moments = mom0;
    for (int d = 0; d < D; ++d) {
      moments[d] * inv_geo_fac;
    }
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Cylindrical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = CoordinateSystems::twopi * mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    auto moments = mom0;
    moments[0] *= CoordinateSystem::twopi *
      (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    // Apply geometry factor
    moments[0] *= inv_geo_fac;
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (axisymmetric) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = CoordinateSystems::twopi * mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    auto moments = mom0;
    moments[0] *= CoordinateSystems::twopi *
      (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[1] *= CoordinateSystems::twopi * rhobar;
    // Apply geometry factor
    for (int d = 0; d < D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (polar) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    auto moments = mom0;
    moments[0] *= (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[1] *= rhobar;
    // Apply geometry factor
    for (int d = 0; d < D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Cylindrical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    mean_radius = 0.5 * (plo[0] + phi[0]);
    auto volume = mean_radius * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // Adjust for different coordinate system
    rhobar = 0.5 * (phi[0] + plo[0]);
    drho_2 = 0.5 * (phi[0] - plo[0]);
    auto moments = mom0;
    moments[0] *= (rhobar*rhobar + drho_2*drho_2/3.0) / rhobar;
    moments[1] *= rhobar;
    moments[2] *= rhobar;
    // Apply geometry factor
    for (int d = 0; d < D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Spherical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
      Point<D> const & plo, Point<D> const & phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (plo[0] + phi[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    auto volume = CoordinateSystems::fourpi *
      (rbar*rbar + dr_2*dr_2/3.0) * vol0;
    // Apply geometry factor
    volume *= inv_geo_fac;
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
    // Adjust for different coordinate system
    rbar = 0.5 * (phi[0] + plo[0]);
    dr_2 = 0.5 * (phi[0] - plo[0]);
    auto moments = mom0;
    moments[0] *= CoordinateSystems::fourpi * (rbar*rbar + dr_2*dr_2);
    // Apply geometry factor
    for (int d = 0; d < D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return moments;
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
  template<int D>
  static constexpr void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Spherical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<int D>
  static constexpr void modify_gradient(Vector<D> & gradient,
      Point<D> const & reference_point) {
    gradient[1] /= reference_point[0];
    gradient[2] /= (reference_point[0] * sin(reference_point[1]));
  }

  /// Modify line element to account for the coordinate system
  template<int D>
  static constexpr void modify_line_element(Vector<D> & line_element,
      Point<D> const & reference_point) {
    line_element[1] *= reference_point[0];
    line_element[2] *= (reference_point[0] * sin(reference_point[1]));
  }

  /// Modify volume to account for the coordinate system
  template<long D>
  static constexpr double modify_volume(double const vol0,
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
    return volume;
  }

  /// Modify moments to account for the coordinate system
  template<long D>
  static constexpr Point<D> modify_moments(Point<D> const & mom0,
      Point<D> const & plo, Point<D> const phi) {
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
    auto moments = mom0;
    moments[0] *= rr1 * ss1;
    moments[1] *= rr2 * (ss1 + cosc_tb * (sinc_dt - cos_dt));
    moments[2] *= rr2 * ss1;
    // Apply geometry factor
    for (int d = 0; d < D; ++d) {
      moments[d] *= inv_geo_fac;
    }
    return moments;
  }

};  // Spherical (3D) Coordinates


// ============================================================================

// Default coordinate system for consistency across the code
using DefaultCoordSys = CartesianCoordinates;

// ============================================================================

}  // namespace Wonton

#endif  // WONTON_SUPPORT_COORDINATESYSTEM_H_
