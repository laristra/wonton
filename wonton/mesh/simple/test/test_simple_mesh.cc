/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/simple/simple_mesh.h"

#include <iostream>
#include <vector>
#include <cmath>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Matrix.h"

#include "gtest/gtest.h"

/*!
  @file test_simple_mesh.cc
  @brief Tests for the simple mesh definition in simple_mesh.h
*/
TEST(Simple_Mesh3D, SingleCell) {
  // Create a single cell mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            1, 1, 1);
  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 3);
  // Cells
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_OWNED),
            1);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);
  // Faces
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            6);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);
  // nodes
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            8);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);

  // These are global id's - the global ordering goes x first, then y, then z
  std::vector<int> nodes;
  // For a cell, the ordering is counterclockwise starting at (x0,y0,z0)
  mesh.cell_get_nodes(0, &nodes);
  // So, `nodes` should have the following global id order:
  std::vector<int> expectedNodes = {0, 1, 3, 2, 4, 5, 7, 6};
  for (int i(0); i < 8; ++i)
    ASSERT_EQ(nodes[i], expectedNodes[i]);

  // Coordinate locations of cell node 5: (1.0, 0.0, 1.0)
  Wonton::Point<3> pp;
  mesh.node_get_coordinates(nodes[5], &pp);
  ASSERT_EQ(1.0, pp[0]);
  ASSERT_EQ(0.0, pp[1]);
  ASSERT_EQ(1.0, pp[2]);
}

TEST(Simple_Mesh3D, SmallGrid) {
  // Create a simple 2x2x2 mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0,
                            1.0, 1.0, 1.0,
                            2, 2, 2);
  // Number of cells
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_OWNED),
            8);
  // Number of faces - interior faces only counted once
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            36);
  // Number of nodes
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            27);

  // Check that the global ID of a node in the center - shared by cells -
  // is the same.
  // This is the global ID of the node in the center of the mesh
  int gid_nodeWant = 13;
  std::vector<int> cellnodes;
  // This is the cell-local index which should correspond to the global index
  std::vector<int> cell_local_index = {6, 7, 5, 4, 2, 3, 1, 0};
  for (int i(0); i < 8; ++i) {
    // Get the nodes associated with this cell
    mesh.cell_get_nodes(i, &cellnodes);
    // Make sure that we match the expected global id
    ASSERT_EQ(gid_nodeWant, cellnodes[cell_local_index[i]]);
  }
}
TEST(Simple_Mesh2D, SingleCell) {
  // Create a single cell mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0,
                            1.0, 1.0,
                              1,   1);
  // Dimensionality
  ASSERT_EQ(mesh.space_dimension(), 2);
  // Cells
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_OWNED),
            1);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);
  // Faces
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            4);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);
  // nodes
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            4);
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_GHOST),
            0);

  // These are global id's - the global ordering goes x first, then y
  std::vector<int> nodes;
  // For a cell, the ordering is counterclockwise starting at (x0,y0)
  mesh.cell_get_nodes(0, &nodes);
  // So, `nodes` should have the following global id order:
  std::vector<int> expectedNodes = {0, 1, 3, 2};
  for (int i(0); i < 4; ++i)
    ASSERT_EQ(nodes[i], expectedNodes[i]);

  // Coordinate locations of cell node 2: (1.0, 1.0)
  Wonton::Point<2> pp;
  mesh.node_get_coordinates(nodes[2], &pp);
  ASSERT_EQ(1.0, pp[0]);
  ASSERT_EQ(1.0, pp[1]);
}

TEST(Simple_Mesh2D, SmallGrid) {
  // Create a simple 2x2x2 mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0,
                            1.0, 1.0,
                              2,   2);
  // Number of cells
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::CELL,
                              Wonton::Entity_type::PARALLEL_OWNED),
            4);
  // Number of faces - interior faces only counted once
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::FACE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            12);
  // Number of nodes
  ASSERT_EQ(mesh.num_entities(Wonton::Entity_kind::NODE,
                              Wonton::Entity_type::PARALLEL_OWNED),
            9);

  // Check that the global ID of a node in the center - shared by cells -
  // is the same.
  // This is the global ID of the node in the center of the mesh
  int gid_nodeWant = 4;
  std::vector<int> cellnodes;
  // This is the cell-local index which should correspond to the global index
  // of the center node for each cell in the mesh
  std::vector<int> cell_local_index = {2,3,1,0};
  for (int i(0); i < 4; ++i) {
    // Get the nodes associated with this cell
    mesh.cell_get_nodes(i, &cellnodes);
    // Make sure that we match the expected global id
    ASSERT_EQ(gid_nodeWant, cellnodes[cell_local_index[i]]);
  }
}


/*!
  @file test_simple_mesh.cc
  @brief Tests for the transform of a simple mesh.
*/
TEST(Simple_Mesh3D, Transform) {
  // Create a single cell mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);

  // rotate Pi/2 radians about the axis {1,1,1} and translate by {10,20,30}
  // by Mathematica
  const double a=1./sqrt(3.), b=1./3.;
  std::vector<std::vector<double>> affinev(3,std::vector<double>(4));
  affinev[0] = {b,b-a,b+a, 10.};
  affinev[1] = {b+a,b,b-a, 20.};
  affinev[2] = {b-a,b+a,b, 30.};
  Wonton::Matrix affine(affinev);

  // do the transform
  mesh.transform<3>(affine);

  // should get these results for the corners
  std::vector<std::vector<double>> expected(8,std::vector<double>(3));
  /*
  expected[0] = {10, 20, 30};
  expected[1] = {10 + b, 20 + a + b, 30 - a + b};
  expected[2] = {10 - a + 2*b, 20 + a + 2*b, 30 + 2*b};
  expected[3] = {10 - a + b, 20 + b, 30 + a + b};
  expected[4] = {10 + a + b, 20 - a + b, 30 + b};
  expected[5] = {10 + a + 2*b, 20 + 2*b, 30 - a + 2*b};
  expected[6] = {10 + 3*b, 20 + 3*b, 30 + 3*b};
  expected[7] = {10 + 2*b, 20 - a + 2*b, 30 + a + 2*b};
  */
  expected[0] = {10, 20, 30};
  expected[1] = {10 + b, 20 + a + b, 30 - a + b};
  expected[2] = {10 - a + b, 20 + b, 30 + a + b};
  expected[3] = {10 - a + 2*b, 20 + a + 2*b, 30 + 2*b};
  expected[4] = {10 + a + b, 20 - a + b, 30 + b};
  expected[5] = {10 + a + 2*b, 20 + 2*b, 30 - a + 2*b};
  expected[6] = {10 + 2*b, 20 - a + 2*b, 30 + a + 2*b};
  expected[7] = {10 + 3*b, 20 + 3*b, 30 + 3*b};

  for (int i=0; i<8; i++) {
    Wonton::Point<3> node;
    mesh.node_get_coordinates(i,&node);
    for (int j=0; j<3; j++) {
      ASSERT_NEAR(node[j], expected[i][j], 1.e-12);
    }
  }
}


/*!
  @file test_simple_mesh.cc
  @brief Tests for the transform of a simple mesh.
*/
TEST(Simple_Mesh2D, Transform) {
  // Create a single cell mesh
  Wonton::Simple_Mesh mesh(0.0, 0.0, 1.0, 1.0, 1, 1);

  // rotate Pi/6 radians and translate by {10,20}
  // by Mathematica
  const double a=sqrt(3.)/2, b=1./2.;
  std::vector<std::vector<double>> affinev(2,std::vector<double>(3));
  affinev[0] = {a,-b, 10.};
  affinev[1] = {b, a, 20.};
  Wonton::Matrix affine(affinev);

  // do the transform
  mesh.transform<2>(affine);

  // should get these results for the corners
  std::vector<std::vector<double>> expected(8,std::vector<double>(3));
  expected[0] = {10, 20};
  expected[1] = {10 + a, 20 + b};
  expected[2] = {10 - b, 20 + a};
  expected[3] = {10 + a - b, 20 + a + b};

  for (int i=0; i<4; i++) {
    Wonton::Point<2> node;
    mesh.node_get_coordinates(i,&node);
    for (int j=0; j<2; j++) {
      ASSERT_NEAR(node[j], expected[i][j], 1.e-12);
    }
  }
}
