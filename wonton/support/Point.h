/*
 This file is part of the Ristra wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

///////////////////////////////////////////////////////////////////////////////
// Copyright 2016 Los Alamos National Laboratory                             //
//                                                                           //
// Original Author: Paul Henning                                             //
//                  TSM, ASCI Problem Setup Team                             //
//                  Applied Physics Division                                 //
//                  Los Alamos National Laboratory                           //
//                  phenning@lanl.gov                                        //
//                                                                           //
// Modified for GK by: Brian Jean                                            //
//                     TSM, ASCI Problem Setup Team                          //
//                     Applied Physics Division                              //
//                     Los Alamos National Laboratory                        //
//                     505.665.6374                                          //
//                     baj@lanl.gov                                          //
//                                                                           //
//                                                                           //
// Modified for Wonton by: Rao Garimella, rao@lanl.gov                       //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef WONTON_POINT_H_
#define WONTON_POINT_H_

#include <cassert>
#include <iostream>
#include <vector>
#include "wonton/support/Vector.h"
#include "wonton/support/wonton.h"

namespace Wonton {

const int X = 0;
const int Y = 1;
const int Z = 2;

/*!
  @class Point "Point.h"
  @brief Represents a point in an N-dimensional space.
  @tparam D Indicates the dimensionality of the Point (this will generally be one
  of [1,2,3]).
*/
template <int D>
class Point {
 private:
  double m_loc[D] {};

 public:

  /// Default constructor - Point at origin in D-space.
  WONTON_INLINE
  Point() {
    for (int i = 0; i < D; i++)
      m_loc[i] = 0.0;
  }

  /*!
    @brief Specialized constructor from a std::vector of arbitary size.
    @param[in] v std::vector of coordinates.
  */
  WONTON_INLINE
  Point(const std::vector<double> &v) {
    assert(v.size() == D);
    for (int i = 0; i < D; i++)
      m_loc[i] = v[i];
  }

  /*!
    @brief Specialized constructor for Points in 3d.
    @param[in] x,y,z The (x,y,z) coordinates of this Point.
  */
  WONTON_INLINE
  Point(const double& x, const double& y, const double& z) {
    assert(D == 3);
    m_loc[0] = x;
    m_loc[1] = y;
    m_loc[2] = z;
  }

  /*!
    @brief Specialized constructor for Points in 2d.
    @param[in] x,y The (x,y) coordinates of this Point.
  */
  WONTON_INLINE
  Point(const double& x, const double& y) {
    assert(D == 2);
    m_loc[0] = x;
    m_loc[1] = y;
  }

  /*!
    @brief Specialized constructor for Points in 1d.
    @param[in] x The x coordinates of this Point.
  */
  WONTON_INLINE
  Point(const double& x) {
    assert(D == 1);
    m_loc[0] = x;
  }

  /// Convert a Vector to a Point.
  WONTON_INLINE
  explicit Point(const Vector<D>& v) {
    for (int i = 0; i < D; i++)
      m_loc[i] = v[i];
  }

  /// Copy constructor.
  WONTON_INLINE
  Point(const Point<D>& rhs) {
    for (int i = 0; i < D; i++)
      m_loc[i] = rhs[i];
  }

  /// Return component @c i of the Point.
  WONTON_INLINE
  const double& operator[](const int& i) const {
    return m_loc[i];
  }

  /// Return component @c i of the Point.
  WONTON_INLINE
  double& operator[](const int& i) {
    return m_loc[i];
  }

  /// Negative of this Point.
  WONTON_INLINE
  Point<D> operator-() const {
    Point<D> p;
    for (int i = 0; i < D; i++) p.m_loc[i] = -m_loc[i];
    return p;
  }

  /// Translate this Point along the Vector @c v.
  WONTON_INLINE
  Point<D>& operator+=(const Vector<D>& v) {
    for (int i = 0; i < D; i++) m_loc[i] += v[i];
    return *this;
  }

  /// Add two Point vectors to get a third point vector
  WONTON_INLINE
  Point<D>& operator+=(const Point<D>& p) {
    for (int i = 0; i < D; i++) m_loc[i] += p[i];
    return *this;
  }

  /// Scale this Point (*)
  WONTON_INLINE
  Point<D>& operator*=(double s) {
    for (int i = 0; i < D; i++) m_loc[i] *= s;
    return *this;
  }

  /// Scale this Point (/)
  WONTON_INLINE
  Point<D>& operator/=(double s) {
    for (int i = 0; i < D; i++) m_loc[i] /= s;
    return *this;
  }

  /// Read in the coordinates a Point from an input stream.
  std::istream& readFromStream(std::istream& is) {
    for (int i = 0; i < D; i++)
      is >> m_loc[i];
    return is;
  }

  /// Pretty printing of the coordinates of a Point to an output stream.
  std::ostream& writeToStream(std::ostream& os) const {
    for (int i = 0; i < D; i++) {
      if (i > 0) os << ' ';
      os << m_loc[i];
    }
    return os;
  }

  /// Convert Point to Vector from the origin to the coordinates of the Point.
  WONTON_INLINE
  Vector<D> asV() {
    Vector<D> v;
    for (int i = 0; i < D; i++)
      v[i] = m_loc[i];
    return v;
  }

  /// Convert Point to Vector from the origin to the coordinates of the Point.
  WONTON_INLINE
  Vector<D> asV() const {
    Vector<D> v;
    for (int i = 0; i < D; i++)
      v[i] = m_loc[i];
    return v;
  }
};

/// Alias for creating a Point in 3d.
typedef Point<3> Point3;
/// Alias for creating a Point in 2d.
typedef Point<2> Point2;
/// Alias for creating a Point in 1d.
typedef Point<1> Point1;

template <int D>
std::ostream& operator<<(std::ostream& os, const Point<D>& p) {
  return p.writeToStream(os);
}


template <int D>
WONTON_INLINE
Point<D> operator+(const Point<D>& p, const Vector<D>& v) {
  return Point<D>(p) += v;
}

template <int D>
WONTON_INLINE
Point<D> operator+(const Point<D>& p1, const Point<D>& p2) {
  return Point<D>(p1) += p2;
}

template <int D>
WONTON_INLINE
Vector<D> operator-(const Point<D>& p1, const Point<D>& p2) {
  Vector<D> v;
  for (int i = 0; i < D; i++) v[i] = p1[i] - p2[i];
  return v;
}

template <int D>
WONTON_INLINE
Point<D> operator*(const Point<D>& p, double s) {
  return Point<D>(p) *= s;
}

template <int D>
WONTON_INLINE
Point<D> operator*(double s, const Point<D>& p) {
  return Point<D>(p) *= s;
}

template <int D>
WONTON_INLINE
Point<D> operator/(const Point<D>& p, double s) {
  return Point<D>(p) /= s;
}

template <int D>
WONTON_INLINE
bool operator==(const Point<D>& p1, const Point<D>& p2) {
  for (int i = 0; i < D; i++)
    if (p1[i] != p2[i])
      return false;
  return true;
}

template <int D>
WONTON_INLINE
bool approxEq(const Point<D>& p1, const Point<D>& p2, double tol = 1.0e-8) {
  // This uses the L_infty norm to avoid nearby subtractions
  for (int i = 0; i < D; i++)
    if (p1[i] < (p2[i] - tol) || p1[i] > (p2[i] + tol))
      return false;
  return true;
}

template <int D>
WONTON_INLINE
bool operator<(const Point<D>& p1, const Point<D>& p2) {
  if (approxEq(p1, p2))
    return false;
  else
    for (int i = 0; i < D; ++i) {
      if (p1[i] < p2[i])
        return true;
      else if (p2[i] < p1[i])
        return false;
    }
  return false;
}


WONTON_INLINE
Point2 ToCylindrical(const Point3& p) {
  Point2 result;

  result[0] = sqrt(p[0]*p[0]+p[1]*p[1]);
  result[1] = p[2];
  return result;
}

WONTON_INLINE
Point3 createP3(double x, double y, double z) {
  Point3 p;

  p[0] = x;
  p[1] = y;
  p[2] = z;

  return p;
}

WONTON_INLINE
Point2 createP2(double x, double y) {
  Point2 p;

  p[0] = x;
  p[1] = y;

  return p;
}

WONTON_INLINE
Point1 createP1(double x) {
  Point1 p;

  p[0] = x;

  return p;
}

}  // namespace Wonton

#endif  // WONTON_POINT_H_
