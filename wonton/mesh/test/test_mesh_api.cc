/*
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/portage/blob/master/LICENSE
*/

#include <iostream>
#include "gtest/gtest.h"

// test includes
#include "my_mesh_wrapper.h"


TEST(Mesh_API, Fixed_D) {
  Wonton::FixedD_Mesh_Wrapper mesh_wrapper;

  // test mesh
  int ncells = mesh_wrapper.num_owned_cells();
  int nfaces = mesh_wrapper.num_owned_faces();
  int nsides = mesh_wrapper.num_owned_sides();
  int nwedges = mesh_wrapper.num_owned_wedges();
  int nnodes = mesh_wrapper.num_owned_nodes();
  int ncorners = mesh_wrapper.num_owned_corners();

  // verify that there is no ghost data
  int ncells_ghost = mesh_wrapper.num_ghost_cells();
  int nfaces_ghost = mesh_wrapper.num_ghost_faces();
  int nnodes_ghost = mesh_wrapper.num_ghost_nodes();
  int isum = ncells_ghost + nfaces_ghost + nnodes_ghost;
  assert(isum == 0);

  assert(ncells == mesh_wrapper.num_entities(Wonton::CELL, Wonton::Entity_type::ALL));
  assert(nfaces == mesh_wrapper.num_entities(Wonton::FACE, Wonton::Entity_type::ALL));
  assert(nsides == mesh_wrapper.num_entities(Wonton::SIDE, Wonton::Entity_type::ALL));
  assert(nnodes == mesh_wrapper.num_entities(Wonton::NODE, Wonton::Entity_type::ALL));

  // cells: accumulate data into arrays and make small calculations to avoid any 
  // possibility for code optimization. 
  std::vector<int> dirs;
  Wonton::Point<3> xc, moment;
  std::vector<std::vector<Wonton::Point<3>>> cell_coords(ncells);
  std::vector<std::vector<int>> cell_cells_f(ncells), cell_cells_n(ncells);
  std::vector<std::vector<int>> cell_sides(ncells), cell_wedges(ncells);
  std::vector<std::vector<int>> cell_faces(ncells), cell_nodes(ncells), cell_corners(ncells);

  std::vector<Wonton::Point<3>> fpoints;
  std::vector<std::array<Wonton::Point<3>, 4>> wpoints;
  std::vector<std::vector<int>> facetpoints;

  for (int c = 0; c < ncells; ++c) {
    mesh_wrapper.cell_get_coordinates(c, &(cell_coords[c]));
    mesh_wrapper.cell_get_faces_and_dirs(c, &(cell_faces[c]), &dirs);
    mesh_wrapper.cell_get_sides(c, &(cell_sides[c]));
    mesh_wrapper.cell_get_wedges(c, &(cell_wedges[c]));
    mesh_wrapper.cell_get_nodes(c, &(cell_nodes[c]));
    mesh_wrapper.cell_get_corners(c, &(cell_corners[c]));
    mesh_wrapper.cell_get_node_adj_cells(c, Wonton::Entity_type::ALL, &(cell_cells_n[c]));
    mesh_wrapper.cell_get_face_adj_cells(c, Wonton::Entity_type::ALL, &(cell_cells_f[c]));

    double vol = mesh_wrapper.cell_volume(c);
    mesh_wrapper.cell_centroid(c, &xc);
    moment += xc * vol;

    mesh_wrapper.cell_get_facetization(c, &facetpoints, &fpoints);
    mesh_wrapper.wedges_get_coordinates(c, &wpoints);
  }

  // faces
  std::vector<std::vector<int>> face_cells(nfaces), face_nodes(nfaces);
  for (int f = 0; f < nfaces; ++f) {
    mesh_wrapper.face_get_cells(f, Wonton::Entity_type::ALL, &(face_cells[f]));
    mesh_wrapper.face_get_nodes(f, &(face_nodes[f]));
  }

  // sides: verify that memory was initialized to zero
  isum = 0;
  double dsum(0.0);
  std::vector<std::array<Wonton::Point<3>, 4>> side_coords(nsides);
  for (int s = 0; s < nsides; ++s) {
    mesh_wrapper.side_get_coordinates(s, &(side_coords[s]));

    int n = mesh_wrapper.side_get_node(s, 0);
    int c = mesh_wrapper.side_get_cell(s);
    int f = mesh_wrapper.side_get_face(s);
    int w = mesh_wrapper.side_get_wedge(s, 0);
    int s2 = mesh_wrapper.side_get_opposite_side(s);
    isum += n + c + f + w + s2;

    double vol = mesh_wrapper.side_volume(s);
    dsum += vol;
  }
  ASSERT_GE(isum, 0);
  ASSERT_LE(fabs(dsum), 1e-20);

  // wedges
  isum = 0;
  dsum = 0.0;
  std::vector<std::array<Wonton::Point<3>, 4>> wedge_coords(nsides);

  for (int w = 0; w < nwedges; ++w) {
    mesh_wrapper.wedge_get_coordinates(w, &(wedge_coords[w]));

    int c = mesh_wrapper.wedge_get_cell(w);
    int s = mesh_wrapper.wedge_get_side(w);
    int cn = mesh_wrapper.wedge_get_corner(w);
    int w2 = mesh_wrapper.wedge_get_opposite_wedge(w);
    int w3 = mesh_wrapper.wedge_get_adjacent_wedge(w);
    int f = mesh_wrapper.wedge_get_face(w);
    int n = mesh_wrapper.wedge_get_node(w);
    isum += c + s + cn + w2 + w3 + f + n;

    double vol = mesh_wrapper.wedge_volume(w);
    dsum += vol;
  }
  ASSERT_GE(isum, 0);
  ASSERT_LE(fabs(dsum), 1e-20);

  // nodes
  std::vector<Wonton::Point<3>> node_coords(nnodes);
  std::vector<std::vector<int>> node_cells(nnodes), node_wedges(nnodes);
  std::vector<std::vector<int>> node_corners(nnodes), node_nodes(nnodes);

  for (int n = 0; n < nnodes; ++n) {
    mesh_wrapper.node_get_coordinates(n, &(node_coords[n]));
    mesh_wrapper.node_get_cells(n, Wonton::Entity_type::ALL, &(node_cells[n]));
    mesh_wrapper.node_get_cell_adj_nodes(n, Wonton::Entity_type::ALL, &(node_nodes[n]));
    mesh_wrapper.node_get_corners(n, Wonton::Entity_type::ALL, &(node_corners[n]));
    mesh_wrapper.node_get_wedges(n, Wonton::Entity_type::ALL, &(node_wedges[n]));
  }

  // corners
  isum = 0;
  dsum = 0.0;
  std::vector<std::vector<int>> corner_wedges(ncorners);

  for (int cn = 0; cn < ncorners; ++cn) {
    mesh_wrapper.corner_get_wedges(cn, &(corner_wedges[cn]));
    int c = mesh_wrapper.corner_get_cell(cn);
    int n = mesh_wrapper.corner_get_node(cn);
    isum += c + n;

    double vol = mesh_wrapper.corner_volume(cn);
    dsum += vol;
  }
  ASSERT_GE(isum, 0);
  ASSERT_LE(fabs(dsum), 1e-20);
}

