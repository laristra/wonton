/*
 This file is part of the Ristra Wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_POLYTOPE_H_
#define WONTON_POLYTOPE_H_

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

namespace Wonton {

template <int D>
class Polytope {
public:
  /*!
    @brief Constructor for a 1D or 2D polygon
    @param vertex_points  Vector of coordinates of the polygon's vertices
    listed counter-clockwise
  */
  Polytope(const std::vector< Point<D> >& vertex_points) {
    assert(D == 1 || D == 2);  // this signature won't work for 3D polytopes
    
    int nfaces = static_cast<int>(vertex_points.size());
    assert((D == 1 && nfaces == 2) || (nfaces >= D+1));

    vertex_points_ = vertex_points;
    face_vertices_.resize(nfaces);
    if (D == 1)
      for (int iface = 0; iface < nfaces; iface++)
        face_vertices_[iface] = { iface };    
    else
      for (int iface = 0; iface < nfaces; iface++)
        face_vertices_[iface] = { iface, (iface + 1)%nfaces };    
  }

  /*!
    @brief Constructor for a general polyhedron
    @param vertex_points  Vector of coordinates of the polyhedron's vertices
    @param face_vertices  Face vertices of the polyhedron

    In 3D, every face is given by IDs of its vertices in
    counter-clockwise order (i.e. so that the normal to the face is
    pointing outside of the polyhedron). In 2D, a "face" is an edge between two 
  */
  Polytope(const std::vector< Point<D> >& vertex_points,
           const std::vector< std::vector<int> >& face_vertices) {
    assert((D == 1 && vertex_points.size() == 2) ||
           (vertex_points.size() >= D+1));
    assert((D == 1 && face_vertices.size() == 2) ||
           (face_vertices.size() >= D+1));
    
    vertex_points_ = vertex_points;
    face_vertices_ = face_vertices;  
  }

  /*! Destructor */
  ~Polytope() = default;

  /*!
    @brief Points for all the vertices of the polytope
    @return  Vector of vertex points
  */
  const std::vector< Point<D> >& points() const { return vertex_points_; }

  /*!
    @brief Indices of vertices of the polytope's face
    in counter-clockwise order (i.e. so that the normal
    to the face is pointing outside of the polytope)
    @param face_id  ID of the face of the polytope
    @return  Vector of indices of the face's vertices
  */
  const std::vector<int>& face_vertices(int const face_id) const {
    assert((face_id >= 0) && (unsigned(face_id) < face_vertices_.size()));
    return face_vertices_[face_id];
  }
  
  /*!
    @brief Coordinates of vertices of the polytope's face
    in counter-clockwise order (i.e. so that the normal
    to the face is pointing outside of the polytope)
    @param face_id  ID of the face of the polytope
    @return  Vector of coordinates of face's vertices
  */
  std::vector< Point<D> > face_points(int const face_id) const {
    assert((face_id >= 0) && (unsigned(face_id) < face_vertices_.size()));

    int nvrts = static_cast<int>(face_vertices_[face_id].size());
    std::vector< Point<D> > fpoints;
    fpoints.reserve(nvrts);
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      fpoints.push_back(vertex_points_[face_vertices_[face_id][ivrt]]);

    return fpoints;
  }

  /*!
    @brief Indices of vertices for all the faces of the polytope
    @return  Vector of vectors of indices of faces' vertices
  */
  const std::vector< std::vector<int> >& face_vertices() const {
    return face_vertices_;
  }
  
  /*!
    @brief Coordinates of a polytope's vertex
    @param vertex_id  ID of a vertex of the polytope
    @return  Coordinates of that vertex
  */
  Point<D> vertex_point(int const vertex_id) const {
    assert((vertex_id >= 0) && (unsigned(vertex_id) < vertex_points_.size()));
    return vertex_points_[vertex_id];
  }
  
  /*!
   @brief Number of material poly's vertices
   @return  Number of vertices
  */
  int num_vertices() const { return static_cast<int>(vertex_points_.size()); }
  /*!
   @brief Number of material poly's faces
   @return  Number of faces
  */
  int num_faces() const { return static_cast<int>(face_vertices_.size()); }

  /*!
    @brief Moments for a face of the polytope
    @param face_id  ID of the face of the polytope
    @return  Vector of moments; moments[0] is the area (length in 2D), moments[i+1]/moments[0]
    is i-th coordinate of the centroid.
  */
  std::vector<double> face_moments(int face_id) const;

  /*!
    @brief Effective normal to a face of the polytope
    @param face_id  ID of the face of the polytope
    @return  Normal vector
  */
  Vector<D> face_normal(int face_id) const;

  /*!
   @brief Moments for the polytope
   @return  Vector of moments; moments[0] is the volume (area in 2D), moments[i+1]/moments[0]
   is i-th coordinate of the centroid
  */  
  std::vector<double> moments(int order = 1) const;

private:
  std::vector< Point<D> > vertex_points_;  // coordinates of vertices
  std::vector< std::vector<int> > face_vertices_;  // vertices of faces
};  // class Polytope

/*!
  @brief Moments for a face of the polygon
  @param face_id  ID of the face of the polygon
  @return  Vector of moments; moments[0] is the size, moments[i+1]/moments[0] is i-th
  coordinate of the centroid
*/
template<>
inline
std::vector<double> Polytope<2>::face_moments(int const face_id) const {
  int nfaces = num_faces();
  assert((face_id >= 0) && (face_id < nfaces));

  std::vector<double> fmoments(3);
  fmoments[0] = (vertex_points_[(face_id + 1)%nfaces] - vertex_points_[face_id]).norm();
  for (int ixy = 0; ixy < 2; ixy++)
    fmoments[ixy + 1] = 0.5*fmoments[0]*(
      vertex_points_[face_id][ixy] + vertex_points_[(face_id + 1)%nfaces][ixy]);

  return fmoments;
}

/*!
  @brief Moments for a face of the polyhedron. Non-triangular faces are assumed to be
  non-flat and are split into triangles using their geometric center.
  @param face_id  ID of the face of the polyhedron
  @return  Vector of moments; moments[0] is the size, moments[i+1]/moments[0] is i-th
  coordinate of the centroid
*/
template<>
inline
std::vector<double> Polytope<3>::face_moments(int const face_id) const {

  assert((face_id >= 0) && (face_id < num_faces()));

  std::vector<double> face_moments(4, 0.0);
  const std::vector<int>& fvertices = face_vertices(face_id);
  int nfvrts = static_cast<int>(fvertices.size());

  if (nfvrts == 3) {
    face_moments[0] = 0.5*cross(vertex_points_[fvertices[1]] - vertex_points_[fvertices[0]], 
                                vertex_points_[fvertices[2]] - vertex_points_[fvertices[0]]).norm();

    for (int ivrt = 0; ivrt < 3; ivrt++)
      for (int ixyz = 0; ixyz < 3; ixyz++)
        face_moments[ixyz + 1] += face_moments[0]*vertex_points_[fvertices[ivrt]][ixyz]/3.0;

    return face_moments;
  }

  Point<3> gcenter;
  for (int ivrt = 0; ivrt < nfvrts; ivrt++)
    gcenter += vertex_points_[fvertices[ivrt]];
  gcenter /= nfvrts;
  
  for (int ivrt = 0; ivrt < nfvrts; ivrt++) {
    int ifv = fvertices[ivrt];
    int isv = fvertices[(ivrt + 1)%nfvrts];

    double tri_size = 0.5*cross(
      vertex_points_[isv] - vertex_points_[ifv], gcenter - vertex_points_[ifv]).norm();

    face_moments[0] += tri_size;
    for (int ixyz = 0; ixyz < 3; ixyz++)
      face_moments[ixyz + 1] += 
        tri_size*(vertex_points_[ifv][ixyz] + vertex_points_[isv][ixyz] + gcenter[ixyz])/3.0;
  }

  return face_moments;
}

/*!
  @brief Normal to a face of the polygon
  @param face_id  ID of the face of the polygon
  @return  Unit normal vector
*/
template<>
inline
Vector<2> Polytope<2>::face_normal(int face_id) const {

  int const nfaces = num_faces();
  assert((face_id >= 0) && (face_id < nfaces));

  int ifv = face_id, isv = (face_id + 1)%nfaces;
  double flen = (vertex_points_[isv] - vertex_points_[ifv]).norm();
  if (flen == 0.0)
    return {0.0, 0.0};

  Vector<2> fnormal(
    (vertex_points_[isv][1] - vertex_points_[ifv][1])/flen,
    -(vertex_points_[isv][0] - vertex_points_[isv][0])/flen);

  return fnormal;
}

/*!
  @brief Effective normal to a face of the polyhedron. Non-triangular 
  faces are assumed to be non-flat and effective normal is computed
  using face centroid for triangulation.
  @param face_id  ID of the face of the polyhedron
  @return  Normal vector, might be non-unitary
*/
template<>
inline
Vector<3> Polytope<3>::face_normal(int face_id) const {

  assert((face_id >= 0) && (face_id < num_faces()));

  Vector<3> fnormal(0.0, 0.0, 0.0);
  const std::vector<int>& fvertices = face_vertices(face_id);
  int nfvrts = static_cast<int>(fvertices.size());
  if (nfvrts == 3) {
    fnormal = cross(vertex_points_[fvertices[1]] - vertex_points_[fvertices[0]], 
                    vertex_points_[fvertices[2]] - vertex_points_[fvertices[0]]);
    fnormal.normalize();
    
    return fnormal;
  }

  std::vector<double> fmoments = face_moments(face_id);
  if (fmoments[0] == 0.0)
    return fnormal;

  Point<3> fcentroid;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    fcentroid[ixyz] = fmoments[ixyz + 1]/fmoments[0];

  for (int ivrt = 0; ivrt < nfvrts; ivrt++) {
    int ifv = fvertices[ivrt], isv = fvertices[(ivrt + 1)%nfvrts];
    Vector<3> tri_normal = cross(
      vertex_points_[isv] - vertex_points_[ifv], 
      fcentroid - vertex_points_[ifv]);
    fnormal += tri_normal;
  }
  fnormal /= fmoments[0];

  return fnormal;
}

/*!
  @brief Moments for the polygon (1D)
  @param moments Computed moments: moments[0] is length, 
  moments[i+1]/moments[0] is i-th coordinate of the centroid, i=1
*/ 
template<>
inline
std::vector<double> Polytope<1>::moments(int order) const {
  assert(order < 4);
  int nmoments = order + 1;
  std::vector<double> poly_moments(nmoments, 0.0);
  poly_moments[0] = vertex_points_[0][1] - vertex_points_[0][0];
  poly_moments[1] = (vertex_points_[0][1] + vertex_points_[0][0]) * poly_moments[0];

  if (order > 1) {
    double b(vertex_points_[0][1]), a(vertex_points_[0][0]);
    double b3(b * b * b), a3(a * a * a);
    poly_moments[2] = (b3 - a3) / 3;
    poly_moments[3] = (b3 * b - a3 * a) / 4;
  }
  return poly_moments;
}

/*!
  @brief Moments for the polygon (2D)
  @param moments Computed moments: moments[0] is area, 
  moments[i+1]/moments[0] is i-th coordinate of the centroid, i=1,2
*/ 
template<>
inline
std::vector<double> Polytope<2>::moments(int order) const {
  assert(order < 3);
  int nmoments = (order + 1) * (order + 2) / 2;
  std::vector<double> poly_moments(nmoments, 0.0);
  int nvrts = num_vertices();

  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    int ifv = ivrt, isv = (ivrt + 1)%nvrts;
    double cur_term = vertex_points_[ifv][0]*vertex_points_[isv][1] - 
                      vertex_points_[ifv][1]*vertex_points_[isv][0];
    poly_moments[0] += cur_term;
    if (order > 0) {
      for (int n = 0; n < 2; ++n)
        poly_moments[n + 1] += cur_term * (vertex_points_[ifv][n] + vertex_points_[isv][n]);
    }
    if (order > 1) {
      double xa = vertex_points_[ifv][0], ya = vertex_points_[ifv][1];
      double xb = vertex_points_[isv][0], yb = vertex_points_[isv][1];
      poly_moments[3] += cur_term * (xa * xa + xa * xb + xb * xb);
      poly_moments[4] += cur_term * (xa * ya + xb * yb + (xa * yb + xb * ya) / 2);
      poly_moments[5] += cur_term * (ya * ya + ya * yb + yb * yb);
    }
  }

  poly_moments[0] /= 2.0;  
  for (int n = 0; n < 2; ++n) poly_moments[n + 1] /= 6.0;
  if (order > 1) {
    for (int n = 0; n < 3; ++n) poly_moments[n + 3] /= 12.0;
  }

  return poly_moments;
}

/*!
  @brief Moments for the polyhedron. Non-triangular faces are assumed to be
  non-flat and are split into triangles using their centroid.
  @param moments Computed moments, moments[0] is volume, 
  moments[i+1]/moments[0] is i-th coordinate of the centroid, i=1,2,3
*/  
template<>
inline
std::vector<double> Polytope<3>::moments(int order) const {
  assert(order < 2);
  std::vector<double> poly_moments(4, 0.0);
  int nfaces = num_faces();

  for (int iface = 0; iface < nfaces; iface++) {
    std::vector<double> fmoments = face_moments(iface);

    if (fmoments[0] == 0.0)
      continue;

    std::vector< Point<3> > face_pts = face_points(iface);
    std::vector< std::vector<int> > itri_pts;
    
    int nfvrts = face_vertices_[iface].size();
    if (nfvrts == 3)
      itri_pts.push_back({0, 1, 2});
    else {
      Point<3> fcentroid;
      for (int ixyz = 0; ixyz < 3; ixyz++)
        fcentroid[ixyz] = fmoments[ixyz + 1]/fmoments[0];
      face_pts.emplace_back(fcentroid);

      itri_pts.reserve(nfvrts);
      for (int ivrt = 0; ivrt < nfvrts; ivrt++)
        itri_pts.push_back({nfvrts, ivrt, (ivrt + 1)%nfvrts});
    }

    for (auto&& point : itri_pts) {
      Vector<3> vcp = cross(face_pts[point[1]] - face_pts[point[0]],
                            face_pts[point[2]] - face_pts[point[0]]);

      poly_moments[0] += dot(vcp, face_pts[point[0]].asV());
      for (int ixyz = 0; ixyz < 3; ixyz++) {
        for (int ivrt = 0; ivrt < 3; ivrt++) {
          int ifv = point[ivrt], isv = point[(ivrt + 1)%3];
          poly_moments[ixyz + 1] += vcp[ixyz]*pow(face_pts[ifv][ixyz] +
                                                  face_pts[isv][ixyz], 2);
        }
      }
    }
  }

  poly_moments[0] /= 6.0;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    poly_moments[ixyz + 1] /= 48.0;

  return poly_moments;
}


}  // namespace Wonton

#endif  // WONTON_POLYTOPE_H_
