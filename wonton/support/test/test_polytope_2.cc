/*
 This file is part of the Ristra Wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/support/Polytope.h"

#include "gtest/gtest.h"
#include "wonton/support/wonton.h"

/// Test the Polytope class for a 2D polygon

TEST(Polytope, Mesh2D) {
  //Test for a unit square
  std::vector< Wonton::Point<2> > square_points = {
    Wonton::Point<2>(0.0, 0.0), Wonton::Point<2>(1.0, 0.0),
    Wonton::Point<2>(1.0, 1.0), Wonton::Point<2>(0.0, 1.0)
  };

  std::vector< std::vector<int> > square_faces = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
  };

  std::vector<double> face_sizes = {1.0, 1.0, 1.0, 1.0};
  std::vector< Wonton::Point<2> > face_centroids = {
    Wonton::Point<2>(0.5, 0.0), Wonton::Point<2>(1.0, 0.5),
    Wonton::Point<2>(0.5, 1.0), Wonton::Point<2>(0.0, 0.5)
  };
  
  //Initialization from ccw ordered vertices
  Wonton::Polytope<2> square_poly(square_points);

  auto const& poly_points    = square_poly.points();
  int const nb_regu_points   = square_points.size();
  int const nb_regu_faces    = square_faces.size();
  int const nb_poly_points   = square_poly.num_vertices();
  int const nb_poly_faces    = square_poly.num_faces();

  ASSERT_EQ(nb_regu_points, nb_poly_points);
  ASSERT_EQ(nb_regu_faces, nb_poly_faces);

  // Verify coordinates
  for (int ivrt = 0; ivrt < nb_regu_points; ivrt++) {
    auto const& regu_point = square_points[ivrt];
    auto const& poly_point = poly_points[ivrt];
    ASSERT_TRUE(approxEq(regu_point, poly_point, 1.0e-15));
  }

  // Verify faces
  for (int iface = 0; iface < nb_regu_faces; iface++) {
    auto const& face_vertices = square_poly.face_vertices(iface);
    ASSERT_EQ(unsigned(2), face_vertices.size());
    ASSERT_EQ(square_faces[iface][0], face_vertices[0]);
    ASSERT_EQ(square_faces[iface][1], face_vertices[1]);
  }
  
  // Verify face moments
  for (int iface = 0; iface < nb_regu_faces; iface++) {
    std::vector<double> regu_moments(3, face_sizes[iface]);
    for (int ixy = 0; ixy < 2; ixy++)
      regu_moments[ixy + 1] *= face_centroids[iface][ixy];
    auto const& poly_moments = square_poly.face_moments(iface);
    for (int im = 0; im < 3; im++)
      ASSERT_NEAR(regu_moments[im], poly_moments[im], 1.0e-15);
  }

  // Test moments for a non-convex pentagon
  std::vector< Wonton::Point<2> > ncv_poly_points = {
    Wonton::Point<2>(0.0, 0.0), Wonton::Point<2>(1.0, 0.0),
    Wonton::Point<2>(1.0, 1.0), Wonton::Point<2>(0.5, 0.5),
    Wonton::Point<2>(0.0, 1.0)
  };

  std::vector<double> ncv_poly_moments = { 0.75, 0.375, 21.0/72.0, 25.0/96.0, 14.0/96.0, 15.0/96.0};
  Wonton::Polytope<2> ncv_poly(ncv_poly_points);

  // Verify moments upto order 2
  auto poly_moments = ncv_poly.moments(2);
  for (int im = 0; im < 6; im++)
    ASSERT_NEAR(ncv_poly_moments[im], poly_moments[im], 1.0e-15);
}
