/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include <numeric>

#include "gtest/gtest.h"

#include "wonton/support/CoordinateSystem.h"

// ============================================================================

constexpr double sqrt2 = 1.414213562373095048;
constexpr double TOLERANCE = 1e-13;

// ============================================================================
// Hard-coded example data

// Get the axis-aligned bounding box for the example cell.
template<typename CoordSys, int D>
auto get_bounding_box();

// Get the reference point for the example calculations.
template<typename CoordSys, int D>
auto get_reference_point();

// Get the gradient before the coordinate system modifies it.
template<typename CoordSys, int D>
auto get_unmodified_gradient();

// Get the gradient after the coordinate system modifies it.
template<typename CoordSys, int D>
auto get_modified_gradient();

// Get the line element before the coordinate system modifies it.
template<typename CoordSys, int D>
auto get_unmodified_line_element();

// Get the line element after the coordinate system modifies it.
template<typename CoordSys, int D>
auto get_modified_line_element();

// Get the volume before it is modified by the coordinate system.
template<typename CoordSys, int D>
auto get_unmodified_volume();

// Get the volume after it is modified by the coordinate system.
template<typename CoordSys, int D>
auto get_modified_volume();

// Get the first moments before they are modified by the coordinate system.
template<typename CoordSys, int D>
auto get_unmodified_first_moments();

// Get the first moments after they are modified by the coordinate system.
template<typename CoordSys, int D>
auto get_modified_first_moments();


// ============================================================================
// Hard-coded example data for 1D Cartesian coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CartesianCoordinates,1>() {
  Wonton::Point<1> lower{1.0};
  Wonton::Point<1> upper{3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CartesianCoordinates,1>() {
  Wonton::Point<1> point{2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CartesianCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CartesianCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CartesianCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CartesianCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CartesianCoordinates,1>() {
  double volume{2.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CartesianCoordinates,1>() {
  double volume{2.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CartesianCoordinates,1>() {
  Wonton::Point<1> first_moments{4.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CartesianCoordinates,1>() {
  Wonton::Point<1> first_moments{4.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for 2D Cartesian coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CartesianCoordinates,2>() {
  Wonton::Point<2> lower{1.0,1.0};
  Wonton::Point<2> upper{3.0,3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CartesianCoordinates,2>() {
  Wonton::Point<2> point{2.0,2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CartesianCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CartesianCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CartesianCoordinates,2>() {
  Wonton::Vector<2> line_element{12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CartesianCoordinates,2>() {
  Wonton::Vector<2> line_element{12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CartesianCoordinates,2>() {
  double volume{4.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CartesianCoordinates,2>() {
  double volume{4.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CartesianCoordinates,2>() {
  Wonton::Point<2> first_moments{8.0,8.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CartesianCoordinates,2>() {
  Wonton::Point<2> first_moments{8.0,8.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for 3D Cartesian coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CartesianCoordinates,3>() {
  Wonton::Point<3> lower{1.0,1.0,1.0};
  Wonton::Point<3> upper{3.0,3.0,3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CartesianCoordinates,3>() {
  Wonton::Point<3> point{2.0,2.0,2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CartesianCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CartesianCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CartesianCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CartesianCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CartesianCoordinates,3>() {
  double volume{8.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CartesianCoordinates,3>() {
  double volume{8.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CartesianCoordinates,3>() {
  Wonton::Point<3> first_moments{16.0,16.0,16.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CartesianCoordinates,3>() {
  Wonton::Point<3> first_moments{16.0,16.0,16.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for cylindrical radial coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Point<1> lower{1.0};
  Wonton::Point<1> upper{3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Point<1> point{2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CylindricalRadialCoordinates,1>() {
  double volume{2.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CylindricalRadialCoordinates,1>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 8.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CylindricalRadialCoordinates,1>() {
  Wonton::Point<1> first_moments{4.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CylindricalRadialCoordinates,1>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<1> first_moments{pi * 52.0 / 3.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for cylindrical axisymmetric coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Point<2> lower{1.0,1.0};
  Wonton::Point<2> upper{3.0,3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Point<2> point{2.0,2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Vector<2> line_element{12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Vector<2> line_element{12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  double volume{4.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 16.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  Wonton::Point<2> first_moments{8.0,8.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CylindricalAxisymmetricCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<2> first_moments{pi * 104.0 / 3.0, pi * 32.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for cylindrical polar coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<2> lower{1.0, pi * 0.5};
  Wonton::Point<2> upper{3.0, pi * 1.5};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<2> point{2.0, pi};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::CylindricalPolarCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::CylindricalPolarCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,6.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::CylindricalPolarCoordinates,2>() {
  Wonton::Vector<2> line_element{12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::CylindricalPolarCoordinates,2>() {
  Wonton::Vector<2> gradient{12.0,24.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 2.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 4.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  Wonton::Point<2> first_moments{pi * 4.0, pi2 * 2.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::CylindricalPolarCoordinates,2>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  Wonton::Point<2> first_moments{pi * 26.0 / 3.0, pi2 * 4.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for cylindrical 3D coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<3> lower{1.0, pi * 0.5, 1.0};
  Wonton::Point<3> upper{3.0, pi * 1.5, 3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<3> point{2.0, pi, 2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::Cylindrical3DCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::Cylindrical3DCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,6.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::Cylindrical3DCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::Cylindrical3DCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,24.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 4.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 8.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  Wonton::Point<3> first_moments{pi * 8.0, pi2 * 4.0, pi * 8.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::Cylindrical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  Wonton::Point<3> first_moments{pi * 52.0 / 3.0, pi2 * 8.0, pi * 16.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for spherical radial coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Point<1> lower{1.0};
  Wonton::Point<1> upper{3.0};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Point<1> point{2.0};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Vector<1> gradient{12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Vector<1> line_element{12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::SphericalRadialCoordinates,1>() {
  double volume{2.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::SphericalRadialCoordinates,1>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 104.0 / 3.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::SphericalRadialCoordinates,1>() {
  Wonton::Point<1> first_moments{4.0};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::SphericalRadialCoordinates,1>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<1> first_moments{pi * 80.0};
  return first_moments;
}


// ============================================================================
// Hard-coded example data for spherical 3D coordinates

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_bounding_box<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<3> lower{1.0, pi * 0.25, pi * 0.5};
  Wonton::Point<3> upper{3.0, pi * 0.75, pi * 1.5};
  return std::make_pair(lower,upper);
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_reference_point<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  Wonton::Point<3> point{2.0, pi * 0.5, pi};
  return point;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_gradient<Wonton::Spherical3DCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,12.0,12.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_gradient<Wonton::Spherical3DCoordinates,3>() {
  Wonton::Vector<3> gradient{12.0,6.0,6.0};
  return gradient;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_line_element<Wonton::Spherical3DCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,12.0,12.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_line_element<Wonton::Spherical3DCoordinates,3>() {
  Wonton::Vector<3> line_element{12.0,24.0,24.0};
  return line_element;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_volume<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  double volume{pi2};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_volume<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  double volume{pi * 26.0 * sqrt2 / 3.0};
  return volume;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_unmodified_first_moments<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi2 = pi * pi;
  constexpr auto pi3 = pi2 * pi;
  Wonton::Point<3> first_moments{pi2 * 2.0, pi3 * 0.5, pi3};
  return first_moments;
}

//  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
template<>
auto get_modified_first_moments<Wonton::Spherical3DCoordinates,3>() {
  using Wonton::CoordinateSystem::pi;
  constexpr auto pi_sqrt2 = pi * sqrt2;
  constexpr auto x = pi_sqrt2 * pi * 13.0 / 3.0;
  Wonton::Point<3> first_moments{pi_sqrt2 * 20.0, x, 2.0 * x};
  return first_moments;
}


// ============================================================================

template<typename CoordSys, int D>
void run_coord_sys_test() {
  // Test verify_coordinate_system
  CoordSys::template verify_coordinate_system<D>();

  // Get a reference point
  Wonton::Point<D> point = get_reference_point<CoordSys,D>();

  // Test modify_gradient
  Wonton::Vector<D> gradient = get_unmodified_gradient<CoordSys,D>();
  CoordSys::modify_gradient(gradient, point);
  Wonton::Vector<D> mod_gradient = get_modified_gradient<CoordSys,D>();
  for (int d = 0; d < D; ++d) {
    ASSERT_NEAR(mod_gradient[d], gradient[d], TOLERANCE);
  }

  // Test modify_line_element
  Wonton::Vector<D> line_element = get_unmodified_line_element<CoordSys,D>();
  CoordSys::modify_line_element(line_element, point);
  Wonton::Vector<D> mod_line_element = get_modified_line_element<CoordSys,D>();
  for (int d = 0; d < D; ++d) {
    ASSERT_NEAR(mod_line_element[d], line_element[d], TOLERANCE);
  }

  // Get box bounds
  auto box = get_bounding_box<CoordSys,D>();
  Wonton::Point<D> lower = std::get<0>(box);
  Wonton::Point<D> upper = std::get<1>(box);

  // Test modify_volume
  double volume = get_unmodified_volume<CoordSys,D>();
  CoordSys::modify_volume(volume, lower, upper);
  double mod_volume = get_modified_volume<CoordSys,D>();
  ASSERT_NEAR(mod_volume, volume, TOLERANCE);

  // Test modify_first_moments
  Wonton::Point<D> first_moments = get_unmodified_first_moments<CoordSys,D>();
  CoordSys::modify_first_moments(first_moments, lower, upper);
  Wonton::Point<D> mod_first_moments = get_modified_first_moments<CoordSys,D>();
  for (int d = 0; d < D; ++d) {
    ASSERT_NEAR(mod_first_moments[d], first_moments[d], TOLERANCE);
  }

  std::cout << CoordSys::str() << std::endl;
}

// ============================================================================

template<int D, int S>
std::vector<double> get_sample_moments();

// ----------------------------------------------------------------------------
// 1D

// Shift by zero (original list)
template<>
std::vector<double> get_sample_moments<1,0>() {
  /* Moments
   * [0] 0, {0}
   * [1] 1, {1}
   * [2] 2, {2}
   * [3] 3, {3}
   */
  return std::vector<double>{0,11,22,33};
}

// Shift by one
template<>
std::vector<double> get_sample_moments<1,1>() {
  /* Moments
   * [0] 0, {0} --> throw away
   * [1] 1, {1} --> 0, {0} [0]
   * [2] 2, {2} --> 1, {1} [1]
   * [3] 3, {3} --> 2, {2} [2]
   */
  return std::vector<double>{11,22,33};
}

// Shift by two
template<>
std::vector<double> get_sample_moments<1,2>() {
  /* Moments
   * [0] 0, {0} --> throw away
   * [1] 1, {1} --> throw away
   * [2] 2, {2} --> 0, {0} [0]
   * [3] 3, {3} --> 1, {1} [1]
   */
  return std::vector<double>{22,33};
}

// ----------------------------------------------------------------------------
// 2D

// Shift by zero (original list)
template<>
std::vector<double> get_sample_moments<2,0>() {
  /* Moments
   * [0] 0, {0,0}
   * [1] 1, {1,0}
   * [2] 1, {0,1}
   * [3] 2, {2,0}
   * [4] 2, {1,1}
   * [5] 2, {0,2}
   * [6] 3, {3,0}
   * [7] 3, {2,1}
   * [8] 3, {1,2}
   * [9] 3, {0,3}
   */
  return std::vector<double>{0,110,101,220,211,202,330,321,312,303};
}

// Shift by one
template<>
std::vector<double> get_sample_moments<2,1>() {
  /* Moments
   * [0] 0, {0,0} --> throw away
   * [1] 1, {1,0} --> 0, {0,0} [0]
   * [2] 1, {0,1} --> throw away
   * [3] 2, {2,0} --> 1, {1,0} [1]
   * [4] 2, {1,1} --> 1, {0,1} [2]
   * [5] 2, {0,2} --> throw away
   * [6] 3, {3,0} --> 2, {2,0} [3]
   * [7] 3, {2,1} --> 2, {1,1} [4]
   * [8] 3, {1,2} --> 2, {0,2} [5]
   * [9] 3, {0,3} --> throw away
   */
  return std::vector<double>{110,220,211,330,321,312};
}

// Shift by two
template<>
std::vector<double> get_sample_moments<2,2>() {
  /* Moments
   * [0] 0, {0,0} --> throw away
   * [1] 1, {1,0} --> throw away
   * [2] 1, {0,1} --> throw away
   * [3] 2, {2,0} --> 0, {0,0} [0]
   * [4] 2, {1,1} --> throw away
   * [5] 2, {0,2} --> throw away
   * [6] 3, {3,0} --> 1, {1,0} [1]
   * [7] 3, {2,1} --> 1, {0,1} [2]
   * [8] 3, {1,2} --> throw away
   * [9] 3, {0,3} --> throw away
   */
  return std::vector<double>{220,330,321};
}

// ----------------------------------------------------------------------------
// 3D

// Shift by zero (original list)
template<>
std::vector<double> get_sample_moments<3,0>() {
  /* Moments
   * [ 0] 0, {0,0,0}
   * [ 1] 1, {1,0,0}
   * [ 2] 1, {0,1,0}
   * [ 3] 1, {0,0,1}
   * [ 4] 2, {2,0,0}
   * [ 5] 2, {1,1,0}
   * [ 6] 2, {1,0,1}
   * [ 7] 2, {0,2,0}
   * [ 8] 2, {0,1,1}
   * [ 9] 2, {0,0,2}
   * [10] 3, {3,0,0}
   * [11] 3, {2,1,0}
   * [12] 3, {2,0,1}
   * [13] 3, {1,2,0}
   * [14] 3, {1,1,1}
   * [15] 3, {1,0,2}
   * [16] 3, {0,3,0}
   * [17] 3, {0,2,1}
   * [18] 3, {0,1,2}
   * [19] 3, {0,0,3}
   */
  return std::vector<double>{0,
    1100,1010,1001,
    2200,2110,2101,2020,2011,2002,
    3300,3210,3201,3120,3111,3102,3030,3021,3012,3003};
}

// Shift by one
template<>
std::vector<double> get_sample_moments<3,1>() {
  /* Moments
   * [ 0] 0, {0,0,0} --> throw away
   * [ 1] 1, {1,0,0} --> 0, {0,0,0} [0]
   * [ 2] 1, {0,1,0} --> throw away
   * [ 3] 1, {0,0,1} --> throw away
   * [ 4] 2, {2,0,0} --> 1, {1,0,0} [1]
   * [ 5] 2, {1,1,0} --> 1, {0,1,0} [2]
   * [ 6] 2, {1,0,1} --> 1, {0,0,1} [3]
   * [ 7] 2, {0,2,0} --> throw away
   * [ 8] 2, {0,1,1} --> throw away
   * [ 9] 2, {0,0,2} --> throw away
   * [10] 3, {3,0,0} --> 2, {2,0,0} [4]
   * [11] 3, {2,1,0} --> 2, {1,1,0} [5]
   * [12] 3, {2,0,1} --> 2, {1,0,1} [6]
   * [13] 3, {1,2,0} --> 2, {0,2,0} [7]
   * [14] 3, {1,1,1} --> 2, {0,1,1} [8]
   * [15] 3, {1,0,2} --> 2, {0,0,2} [9]
   * [16] 3, {0,3,0} --> throw away
   * [17] 3, {0,2,1} --> throw away
   * [18] 3, {0,1,2} --> throw away
   * [19] 3, {0,0,3} --> throw away
   */
  return
    std::vector<double>{1100,2200,2110,2101,3300,3210,3201,3120,3111,3102};
}

// Shift by two
template<>
std::vector<double> get_sample_moments<3,2>() {
  /* Moments
   * [ 0] 0, {0,0,0} --> throw away
   * [ 1] 1, {1,0,0} --> throw away
   * [ 2] 1, {0,1,0} --> throw away
   * [ 3] 1, {0,0,1} --> throw away
   * [ 4] 2, {2,0,0} --> 0, {0,0,0} [0]
   * [ 5] 2, {1,1,0} --> throw away
   * [ 6] 2, {1,0,1} --> throw away
   * [ 7] 2, {0,2,0} --> throw away
   * [ 8] 2, {0,1,1} --> throw away
   * [ 9] 2, {0,0,2} --> throw away
   * [10] 3, {3,0,0} --> 1, {1,0,0} [1]
   * [11] 3, {2,1,0} --> 1, {0,1,0} [2]
   * [12] 3, {2,0,1} --> 1, {0,0,1} [3]
   * [13] 3, {1,2,0} --> throw away
   * [14] 3, {1,1,1} --> throw away
   * [15] 3, {1,0,2} --> throw away
   * [16] 3, {0,3,0} --> throw away
   * [17] 3, {0,2,1} --> throw away
   * [18] 3, {0,1,2} --> throw away
   * [19] 3, {0,0,3} --> throw away
   */
  return
    std::vector<double>{2200,3300,3210,3201};
}

// ============================================================================
// This needs to be separated out from run_coord_sys_test because not all
// coordinate systems enable this feature.

template<typename CoordSys, int D>
constexpr void run_shift_test() {
  // Construct original list of moments
  auto moments = get_sample_moments<D,0>();
  // Shift moments
  //    For many coordinate systems, a moment in that coordinate system is
  // equal to a higher-order moment in Cartesian coordinates.  Portage does
  // many of its calculations assuming Cartesian coordinates, and this method
  // will shift the moments around to get the correct moments in the current
  // coordinate system.
  CoordSys::template shift_moments_list<D>(moments);
  // Construct expected list of moments
  auto new_moments = get_sample_moments<D,CoordSys::moment_shift>();
  // Verify
  ASSERT_EQ(new_moments.size(), moments.size());
  for (int d = 0; d < new_moments.size(); ++d) {
    ASSERT_NEAR(new_moments[d],
            moments[d] / CoordSys::moment_coefficient, TOLERANCE);
  }
}

// ============================================================================

// Cartesian 1D
TEST(Moment_Index_test, Cartesian1D) {
  using CoordSys = Wonton::CartesianCoordinates;
  constexpr int D = 1;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// Cartesian 2D
TEST(Moment_Index_test, Cartesian2D) {
  using CoordSys = Wonton::CartesianCoordinates;
  constexpr int D = 2;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// Cartesian 3D
TEST(Moment_Index_test, Cartesian3D) {
  using CoordSys = Wonton::CartesianCoordinates;
  constexpr int D = 3;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// cylindrical radial
TEST(Moment_Index_test, CylindricalRadial) {
  using CoordSys = Wonton::CylindricalRadialCoordinates;
  constexpr int D = 1;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// cylindrical axisymmetric
TEST(Moment_Index_test, CylindricalAxisymmetric) {
  using CoordSys = Wonton::CylindricalAxisymmetricCoordinates;
  constexpr int D = 2;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// cylindrical polar
TEST(Moment_Index_test, CylindricalPolar) {
  using CoordSys = Wonton::CylindricalPolarCoordinates;
  constexpr int D = 2;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// cylindrical 3D
TEST(Moment_Index_test, Cylindrical3D) {
  using CoordSys = Wonton::Cylindrical3DCoordinates;
  constexpr int D = 3;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// spherical radial
TEST(Moment_Index_test, SphericalRadial) {
  using CoordSys = Wonton::SphericalRadialCoordinates;
  constexpr int D = 1;
  run_coord_sys_test<CoordSys,D>();
  run_shift_test<CoordSys,D>();
}

// spherical 3D
TEST(Moment_Index_test, Spherical3D) {
  using CoordSys = Wonton::Spherical3DCoordinates;
  constexpr int D = 3;
  run_coord_sys_test<CoordSys,D>();
  // Spherical 3D coordinates is not permitted to use the moment-shift method.
  //run_shift_test<CoordSys,D>();
}

