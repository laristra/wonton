/*
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/
#include <gtest/gtest.h>
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"

TEST(Vector, BasicConstruct) {

  // default initialization
  Wonton::Vector<2> v0;
  ASSERT_DOUBLE_EQ(v0[0], 0.0);
  ASSERT_DOUBLE_EQ(v0[0], v0[1]);

  // initializing from a standard vector
  std::vector<double> l1 = { 0.4 };
  std::vector<double> l2 = { 0.1, 2.3 };
  std::vector<double> l3 = { 0.1, 2.3, 3.5 };
  Wonton::Vector<1> v1(l1);
  Wonton::Vector<2> v2(l2);
  Wonton::Vector<3> v3(l3);
  ASSERT_DOUBLE_EQ(v1[0], l1[0]);
  ASSERT_DOUBLE_EQ(v2[0], l2[0]);
  ASSERT_DOUBLE_EQ(v2[1], l2[1]);
  ASSERT_DOUBLE_EQ(v3[0], v2[0]);
  ASSERT_DOUBLE_EQ(v3[1], v2[1]);
  ASSERT_DOUBLE_EQ(v3[2], l3[2]);

  // direct initialization
  Wonton::Vector<1> v4(1.2);
  Wonton::Vector<2> v5(1.2, 3.4);
  Wonton::Vector<3> v6(1.2, 3.4, 4.6);
  ASSERT_DOUBLE_EQ(v4[0], 1.2);
  ASSERT_DOUBLE_EQ(v5[0], v4[0]);
  ASSERT_DOUBLE_EQ(v5[1], 3.4);
  ASSERT_DOUBLE_EQ(v6[0], v5[0]);
  ASSERT_DOUBLE_EQ(v6[1], v5[1]);
  ASSERT_DOUBLE_EQ(v6[2], 4.6);

  // copy and assign
  Wonton::Vector<2> v7(v5);
  Wonton::Vector<2> v8 = v7;
  ASSERT_DOUBLE_EQ(v7[0], v5[0]);
  ASSERT_DOUBLE_EQ(v7[1], v5[1]);
  ASSERT_DOUBLE_EQ(v8[0], v7[0]);
  ASSERT_DOUBLE_EQ(v8[1], v7[1]);
}

TEST(Vector, ImplicitConstruct) {

  // copy-list-initialization
  Wonton::Vector<1> v1 = {0.5};
  Wonton::Vector<1> v2 = {0.5, 1.3};
  ASSERT_DOUBLE_EQ(v1[0], 0.5);
  ASSERT_DOUBLE_EQ(v2[0], v1[0]);
  ASSERT_DOUBLE_EQ(v2[1], 1.3);

  // list-initialization
  Wonton::Vector<1> v3 {0.4};
  Wonton::Vector<2> v4 {0.1, 2.3};
  ASSERT_DOUBLE_EQ(v3[0], 0.4);
  ASSERT_DOUBLE_EQ(v4[0], 0.1);
  ASSERT_DOUBLE_EQ(v4[1], 2.3);

  // collection of points from an initializer list
  std::vector<Wonton::Vector<2>> list = {{0.4, 0.5},
                                         {1.3, 2.5}};

  ASSERT_DOUBLE_EQ(list[0][0], 0.4);
  ASSERT_DOUBLE_EQ(list[0][1], 0.5);
  ASSERT_DOUBLE_EQ(list[1][0], 1.3);
  ASSERT_DOUBLE_EQ(list[1][1], 2.5);
}

TEST(Vector, Operators) {

  Wonton::Vector<2> v(1.2, 3.4);
  Wonton::Vector<2> q = v;
  Wonton::Vector<2> r = v;
  Wonton::Vector<2> w = v;
  Wonton::Vector<2> t(v[0] - 1.E-9, v[1]);

  // addition
  ASSERT_DOUBLE_EQ((v + q)[0], 2 * v[0]);
  ASSERT_DOUBLE_EQ((v + q)[1], 2 * v[1]);

  // negation
  ASSERT_DOUBLE_EQ((v - q)[0], 0.0);
  ASSERT_DOUBLE_EQ((v - q)[1], 0.0);
  ASSERT_DOUBLE_EQ(-v[0], -1.2);
  ASSERT_DOUBLE_EQ(-v[1], -3.4);

  // incrementation
  q += v;
  ASSERT_DOUBLE_EQ(q[0], v[0] + r[0]);
  ASSERT_DOUBLE_EQ(q[1], v[1] + r[1]);

  r += q;
  ASSERT_DOUBLE_EQ(r[0], 3 * v[0]);
  ASSERT_DOUBLE_EQ(r[1], 3 * v[1]);

  std::vector<double> l = { v[0], v[1] };
  r = v;
  r += l;
  ASSERT_DOUBLE_EQ(r[0], v[0] + l[0]);
  ASSERT_DOUBLE_EQ(r[1], v[1] + l[1]);

  // scaling and scalar product
  double const s = 10.;
  ASSERT_DOUBLE_EQ((v * s)[0], v[0] * s);
  ASSERT_DOUBLE_EQ((v * s)[1], v[1] * s);

  w *= s;
  ASSERT_DOUBLE_EQ(w[0], v[0] * s);
  ASSERT_DOUBLE_EQ(w[1], v[1] * s);
  ASSERT_DOUBLE_EQ((w / s)[0], v[0]);
  ASSERT_DOUBLE_EQ((w / s)[1], v[1]);

  w /= s;
  ASSERT_DOUBLE_EQ(w[0], v[0]);
  ASSERT_DOUBLE_EQ(w[1], v[1]);

  // L^[1|2|inf] norm
  ASSERT_DOUBLE_EQ(v.one_norm(), std::fabs(v[0]) + std::fabs(v[1]));
  ASSERT_DOUBLE_EQ(v.norm(), sqrt(v[0]*v[0] + v[1]*v[1]));
  ASSERT_DOUBLE_EQ(v.max_norm(), std::max(std::fabs(v[0]), std::fabs(v[1])));

  // normalize and zero
  Wonton::Vector<2> dummy = v;
  dummy.normalize();
  ASSERT_DOUBLE_EQ(dummy.norm(), 1.0);

  // set values
  dummy.zero();
  ASSERT_DOUBLE_EQ(dummy[0], 0.0);
  ASSERT_DOUBLE_EQ(dummy[1], 0.0);
  ASSERT_DOUBLE_EQ(dummy.norm(), 0.0);
  ASSERT_TRUE(dummy.is_zero(1.E-10));

  dummy.fill(2.0);
  ASSERT_DOUBLE_EQ(dummy[0], 2.0);
  ASSERT_DOUBLE_EQ(dummy[1], dummy[0]);

  dummy.axis(1);
  ASSERT_DOUBLE_EQ(dummy[0], 0.0);
  ASSERT_DOUBLE_EQ(dummy[1], 1.0);

  // dot product, determinant and max component
  int i = -1;
  ASSERT_DOUBLE_EQ(Wonton::dot(v, q), v[0] * q[0] + v[1] * q[1]);
  ASSERT_DOUBLE_EQ(Wonton::dot(v, v), std::pow(v.norm(), 2));
  ASSERT_DOUBLE_EQ(Wonton::cross(v, q), v[0] * q[1] - q[0] * v[1]);
  ASSERT_DOUBLE_EQ(Wonton::MaxComponent(v, i), std::max(v[0], v[1]));

  // cross product
  Wonton::Vector<3> x(1.,0.,0.);
  Wonton::Vector<3> y(0.,1.,0.);
  Wonton::Vector<3> z = Wonton::cross(x, y);
  ASSERT_DOUBLE_EQ(std::abs(z[0]), 0.0);
  ASSERT_DOUBLE_EQ(std::abs(z[1]), 0.0);
  ASSERT_DOUBLE_EQ(std::abs(z[2]), 1.0);

}