/*
 This file is part of the Ristra Wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/support/Polytope.h"

#include "gtest/gtest.h"
#include "wonton/support/wonton.h"

/// Test the Polytope class for a 3D polyhedron

TEST(Polytope, Mesh3D) {
  //Test for a right triangular prism
  std::vector< Wonton::Point<3> > prism_points = {
    Wonton::Point<3>(1.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
    Wonton::Point<3>(0.0, 0.0, 0.0), Wonton::Point<3>(0.0, 1.0, 1.0),
    Wonton::Point<3>(1.0, 0.0, 1.0), Wonton::Point<3>(0.0, 0.0, 1.0)
  };

  std::vector< std::vector<int> > prism_faces = {
    {0, 2, 1}, {2, 0, 4, 5}, {3, 4, 0, 1}, {5, 3, 1, 2}, {3, 5, 4}
  };

  std::vector<double> prism_face_sizes = {0.5, 1.0, sqrt(2), 1.0, 0.5};
  std::vector< Wonton::Point<3> > face_centroids = {
    Wonton::Point<3>(1.0/3.0, 1.0/3.0, 0.0),
    Wonton::Point<3>(0.5, 0.0, 0.5),
    Wonton::Point<3>(0.5, 0.5, 0.5),
    Wonton::Point<3>(0.0, 0.5, 0.5),
    Wonton::Point<3>(1.0/3.0, 1.0/3.0, 1.0)
  };
  
  //Initialization
  Wonton::Polytope<3> prism_poly(prism_points, prism_faces);
  
  //Verify coordinates
  auto const& poly_points    = prism_poly.points();
  int const nb_prism_points  = prism_points.size();
  int const nb_prism_faces   = prism_faces.size();
  int const nb_poly_points   = prism_poly.num_vertices();
  int const nb_poly_faces    = prism_poly.num_faces();

  ASSERT_EQ(nb_prism_points, nb_poly_points);
  for (int ivrt = 0; ivrt < nb_prism_points; ivrt++)
    ASSERT_TRUE(approxEq(prism_points[ivrt], poly_points[ivrt], 1.0e-15));
  
  //Verify faces
  ASSERT_EQ(nb_prism_faces, nb_poly_faces);
  for (int iface = 0; iface < nb_prism_faces; iface++) {
    auto const& face_vertices = prism_poly.face_vertices(iface);
    int const nb_regu_face_points = prism_faces[iface].size();
    int const nb_poly_face_points = face_vertices.size();
    ASSERT_EQ(nb_regu_face_points, nb_poly_face_points);

    for (int ivrt = 0; ivrt < nb_regu_face_points; ivrt++)
      ASSERT_EQ(prism_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify face moments
  for (int iface = 0; iface < nb_prism_faces; iface++) {
    std::vector<double> face_moments(4, prism_face_sizes[iface]);
    for (int ixyz = 0; ixyz < 3; ixyz++)
      face_moments[ixyz + 1] *= face_centroids[iface][ixyz];    
    auto const& poly_moments = prism_poly.face_moments(iface);
    for (int im = 0; im < 4; im++)
      ASSERT_NEAR(face_moments[im], poly_moments[im], 1.0e-15);
  }
  
  //Non-convex distorted prism
  std::vector< Wonton::Point<3> > ncv_prism_points = {
    Wonton::Point<3>(1.0, 0.0, 0.0), Wonton::Point<3>(0.4, 0.8, 0.2),
    Wonton::Point<3>(1.0, 1.0, 0.0), Wonton::Point<3>(0.0, 1.0, 0.0),
    Wonton::Point<3>(0.5, 0.1, 0.1), Wonton::Point<3>(0.0, 0.0, 0.0),
    Wonton::Point<3>(0.0, 0.0, 1.0), Wonton::Point<3>(1.0, 1.0, 1.0),
    Wonton::Point<3>(0.5, 0.9, 1.1), Wonton::Point<3>(1.0, 0.0, 1.0),
    Wonton::Point<3>(0.0, 1.0, 1.0), Wonton::Point<3>(0.6, 0.2, 0.8)
  };

  std::vector< std::vector<int> > ncv_prism_faces = {
    {3, 1, 2, 0, 4, 5},
    {5, 4, 11, 6}, {0, 9, 11, 4}, {7, 9, 0, 2},
    {2, 1, 8, 7}, {8, 1, 3, 10}, {5, 6, 10, 3},
    {7, 8, 10, 6, 11, 9}
  };

  std::vector<double> ncv_prism_face_sizes = {
    0.90144108856477323, 
    0.488981965269429497, 0.406051769001555596, 1.0, 
    0.551936491486091696, 0.449899641235988745, 1.0, 
    0.917588789851955799
  };
  std::vector<Wonton::Point<3>> ncv_prism_face_centroids = {
    Wonton::Point<3>(0.49982029799691324, 0.48793283780407681, 0.038845671672909921),
    Wonton::Point<3>(0.26077565588523149, 0.071880513370773835, 0.49415030095918006),
    Wonton::Point<3>(0.78729807432173626, 0.070061777159007799, 0.46574083274162237),
    Wonton::Point<3>(1.0, 0.5, 0.5),
    Wonton::Point<3>(0.72714164844017881, 0.92490808257951729, 0.55650182505587276),
    Wonton::Point<3>(0.22120243710721832, 0.92700420108481663, 0.58407100072352491),
    Wonton::Point<3>(0.0, 0.5, 0.5),
    Wonton::Point<3>(0.50063253857544732, 0.51180586493539493, 0.98662997695018306)
  };

  std::vector<double> ncv_prism_moments = {
    0.80774747387852563, 0.40401879823038706,
    0.40692862723931095, 0.41298799905686806
  };
  
  //Initialization
  Wonton::Polytope<3> ncv_prism_poly(ncv_prism_points, ncv_prism_faces);

  int const nb_ncv_prism_faces = ncv_prism_faces.size();

  //Verify face moments
  for (int iface = 0; iface < nb_ncv_prism_faces; iface++) {
    std::vector<double> face_moments(4, ncv_prism_face_sizes[iface]);
    for (int ixyz = 0; ixyz < 3; ixyz++)
      face_moments[ixyz + 1] *= ncv_prism_face_centroids[iface][ixyz];    
    auto const& poly_moments = ncv_prism_poly.face_moments(iface);
    for (int im = 0; im < 4; im++)
      ASSERT_NEAR(face_moments[im], poly_moments[im], 1.0e-15);
  }

  //Verify moments
  auto matpoly_moments = ncv_prism_poly.moments();
  for (int im = 0; im < 4; im++)
    ASSERT_NEAR(ncv_prism_moments[im], matpoly_moments[im], 1.0e-15);    
}
