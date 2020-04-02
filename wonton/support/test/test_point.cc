/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/
#include <gtest/gtest.h>
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"


TEST(Point, BasicConstruct) {

  // default initialization
  Wonton::Point<2> p0;
  ASSERT_DOUBLE_EQ(p0[0], 0.0);
  ASSERT_DOUBLE_EQ(p0[0], p0[1]);

  // initializing from a standard vector
  std::vector<double> v1 = { 0.4 };
  std::vector<double> v2 = { 0.1, 2.3 };
  std::vector<double> v3 = { 0.1, 2.3, 3.5 };
  Wonton::Point<1> p1(v1);
  Wonton::Point<2> p2(v2);
  Wonton::Point<3> p3(v3);
  ASSERT_DOUBLE_EQ(p1[0], v1[0]);
  ASSERT_DOUBLE_EQ(p2[0], v2[0]);
  ASSERT_DOUBLE_EQ(p2[1], v2[1]);
  ASSERT_DOUBLE_EQ(p3[0], p2[0]);
  ASSERT_DOUBLE_EQ(p3[1], p2[1]);
  ASSERT_DOUBLE_EQ(p3[2], v3[2]);

  // direct initialization
  Wonton::Point<1> p4(1.2);
  Wonton::Point<2> p5(1.2, 3.4);
  Wonton::Point<3> p6(1.2, 3.4, 4.6);
  ASSERT_DOUBLE_EQ(p4[0], 1.2);
  ASSERT_DOUBLE_EQ(p5[0], p4[0]);
  ASSERT_DOUBLE_EQ(p5[1], 3.4);
  ASSERT_DOUBLE_EQ(p6[0], p5[0]);
  ASSERT_DOUBLE_EQ(p6[1], p5[1]);
  ASSERT_DOUBLE_EQ(p6[2], 4.6);

  // copy and assign
  Wonton::Point<2> p7(p5);
  Wonton::Point<2> p8 = p7;
  ASSERT_DOUBLE_EQ(p7[0], p5[0]);
  ASSERT_DOUBLE_EQ(p7[1], p5[1]);
  ASSERT_DOUBLE_EQ(p8[0], p7[0]);
  ASSERT_DOUBLE_EQ(p8[1], p7[1]);

  // use dedicated methods
  Wonton::Point<1> q1 = Wonton::createP1(p1[0]);
  Wonton::Point<2> q2 = Wonton::createP2(p2[0], p2[1]);
  Wonton::Point<3> q3 = Wonton::createP3(p3[0], p3[1], p3[2]);
  ASSERT_DOUBLE_EQ(q1[0], p1[0]);
  ASSERT_DOUBLE_EQ(q2[0], p2[0]);
  ASSERT_DOUBLE_EQ(q2[1], p2[1]);
  ASSERT_DOUBLE_EQ(q3[0], p3[0]);
  ASSERT_DOUBLE_EQ(q3[1], p3[1]);
  ASSERT_DOUBLE_EQ(q3[2], p3[2]);
}

TEST(Point, ImplicitConstruct) {

  // copy-list-initialization
  Wonton::Point<1> p1 = {0.5};
  Wonton::Point<1> p2 = {0.5, 1.3};
  ASSERT_DOUBLE_EQ(p1[0], 0.5);
  ASSERT_DOUBLE_EQ(p2[0], p1[0]);
  ASSERT_DOUBLE_EQ(p2[1], 1.3);

  // list-initialization
  Wonton::Point<1> p3 {0.4};
  Wonton::Point<2> p4 {0.1, 2.3};
  ASSERT_DOUBLE_EQ(p3[0], 0.4);
  ASSERT_DOUBLE_EQ(p4[0], 0.1);
  ASSERT_DOUBLE_EQ(p4[1], 2.3);

  // collection of points from an initializer list
  std::vector<Wonton::Point<2>> points = {{0.4, 0.5},
                                          {1.3, 2.5}};

  ASSERT_EQ(points.size(), unsigned(2));
  ASSERT_DOUBLE_EQ(points[0][0], 0.4);
  ASSERT_DOUBLE_EQ(points[0][1], 0.5);
  ASSERT_DOUBLE_EQ(points[1][0], 1.3);
  ASSERT_DOUBLE_EQ(points[1][1], 2.5);
}

TEST(Point, Operators) {

  Wonton::Point<2> p(1.2, 3.4);
  Wonton::Point<2> q = p;
  Wonton::Point<2> r = p;
  Wonton::Point<2> z = p;
  Wonton::Point<2> t(p[0] - 1.E-9, p[1]);
  Wonton::Vector<2> v(p[0], p[1]);

  // equality
  ASSERT_TRUE(p == q);
  ASSERT_TRUE(q == r);
  ASSERT_TRUE(r == z);
  ASSERT_TRUE(Wonton::approxEq(p, q));
  ASSERT_TRUE(Wonton::approxEq(p, t));
  ASSERT_FALSE(Wonton::approxEq(p, t, 1.E-12));

  // inequality
  Wonton::Point<2> a(-1.1, 3.4);
  Wonton::Point<2> b(a[1], a[0]);
  ASSERT_FALSE(p < q);
  ASSERT_TRUE(a < q);
  ASSERT_TRUE(a < b);
  ASSERT_FALSE(b < a);

  // addition
  ASSERT_DOUBLE_EQ((p + q)[0], 2 * p[0]);
  ASSERT_DOUBLE_EQ((p + q)[1], 2 * p[1]);
  ASSERT_DOUBLE_EQ((p + v)[0], 2 * p[0]);
  ASSERT_DOUBLE_EQ((p + v)[1], 2 * p[1]);

  // negation
  ASSERT_DOUBLE_EQ((p - q)[0], 0.0);
  ASSERT_DOUBLE_EQ((p - q)[1], 0.0);
  ASSERT_DOUBLE_EQ(-p[0], -1.2);
  ASSERT_DOUBLE_EQ(-p[1], -3.4);

  // incrementation
  q += v;
  ASSERT_DOUBLE_EQ(q[0], p[0] + v[0]);
  ASSERT_DOUBLE_EQ(q[1], p[1] + v[1]);

  r += q;
  ASSERT_DOUBLE_EQ(r[0], 2 * p[0] + v[0]);
  ASSERT_DOUBLE_EQ(r[1], 2 * p[1] + v[1]);

  std::vector<double> l = { p[0], p[1] };
  r = p;
  r += l;
  ASSERT_DOUBLE_EQ(r[0], p[0] + l[0]);
  ASSERT_DOUBLE_EQ(r[1], p[1] + l[1]);

  // scaling
  double const s = 10.;
  ASSERT_DOUBLE_EQ((p * s)[0], p[0] * s);
  ASSERT_DOUBLE_EQ((p * s)[1], p[1] * s);

  z *= s;
  ASSERT_DOUBLE_EQ(z[0], p[0] * s);
  ASSERT_DOUBLE_EQ(z[1], p[1] * s);
  ASSERT_DOUBLE_EQ((z / s)[0], p[0]);
  ASSERT_DOUBLE_EQ((z / s)[1], p[1]);

  z /= s;
  ASSERT_DOUBLE_EQ(z[0], p[0]);
  ASSERT_DOUBLE_EQ(z[1], p[1]);

  // conversion
  Wonton::Vector<2> asv = p.asV();
  ASSERT_DOUBLE_EQ(asv[0], p[0]);
  ASSERT_DOUBLE_EQ(asv[1], p[1]);
}