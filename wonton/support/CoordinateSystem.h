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
 */

namespace Wonton {

// ============================================================================
/// Cartesian Coordinates
struct CartesianCoordinates {
 public:

  /// Geometry factor
  static constexpr double geometry_factor = 1;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid for any positive dimensionality
    static_assert(D >= 1,
        "Cartesian coordinates must have positive dimensionality.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(line_element));
  }
};  // Cartesian Coordinates


// ============================================================================
/// Cylindrical (Radial) Coordinates
/// Only valid in 1D.  The coordinate is the distance from the z axis.
struct CylindricalRadialCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 2 pi.  EAP uses pi.
  static constexpr double geometry_factor = 1;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Cylindrical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(line_element));
  }
};  // Cylindrical (Radial) Coordinates


// ============================================================================
/// Cylindrical (Axisymmetric) Coordinates
/// Only valid in 2D.  The coordinates are the distance from the z axis and the
/// height.
struct CylindricalAxisymmetricCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 2 pi.  EAP uses pi.
  static constexpr double geometry_factor = 1;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (axisymmetric) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(line_element));
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

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 2D
    static_assert(D == 2,
        "Cylindrical (polar) coordinates only valid in 2D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    auto new_gradient = gradient;
    new_gradient[1] /= reference_point[0];
    return(std::move(new_gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    auto new_line_element = line_element;
    new_line_element[1] *= reference_point[0];
    return(std::move(new_line_element));
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

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Cylindrical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    auto new_gradient = gradient;
    new_gradient[1] /= reference_point[0];
    return(std::move(new_gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    auto new_line_element = line_element;
    new_line_element[1] *= reference_point[0];
    return(std::move(new_line_element));
  }
};  // Cylindrical (3D) Coordinates


// ============================================================================
/// Spherical (Radial) Coordinates
/// Only valid in 1D.  The coordinate is the distance from the origin.
struct SphericalRadialCoordinates {
 public:

  /// Geometry factor
  /// A very common choice is 4 pi.  EAP uses 4 pi / 3.
  static constexpr double geometry_factor = 1;

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 1D
    static_assert(D == 1,
        "Spherical (radial) coordinates only valid in 1D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    // No change from "standard", Cartesian-like calculation.
    return(std::move(line_element));
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

  /// Verify coordinate system / dimensionality combination
  template<long D>
  static void verify_coordinate_system() {
    // Valid only in 3D
    static_assert(D == 3,
        "Spherical (3D) coordinates only valid in 3D.");
  }

  /// Modify gradient to account for the coordinate system
  template<long D>
  static Vector<D> modify_gradient(Vector<D> const & gradient,
      Point<D> const & reference_point) {
    auto new_gradient = gradient;
    new_gradient[1] /= reference_point[0];
    new_gradient[2] /= (reference_point[0] * sin(reference_point[1]));
    return(std::move(new_gradient));
  }

  /// Modify line element to account for the coordinate system
  template<long D>
  static Vector<D> modify_line_element(Vector<D> const & line_element,
      Point<D> const & reference_point) {
    auto new_line_element = line_element;
    new_line_element[1] *= reference_point[0];
    new_line_element[2] *= (reference_point[0] * sin(reference_point[1]));
    return(std::move(new_line_element));
  }
};  // Spherical (3D) Coordinates


// ============================================================================

// Default coordinate system for consistency across the code
using DefaultCoordSys = CartesianCoordinates;

// ============================================================================

}  // namespace Wonton

#endif  // WONTON_SUPPORT_COORDINATESYSTEM_H_
