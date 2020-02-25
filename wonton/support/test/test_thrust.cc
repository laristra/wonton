/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include "thrust/transform.h"

TEST(THRUST_TEST, TRANSFORM) {
  Wonton::counting_iterator cellids_begin(0);
  Wonton::counting_iterator cellids_end(4);
  ASSERT_EQ(typeid(cellids_begin),
            typeid(thrust::counting_iterator<unsigned int>));
  ASSERT_EQ(typeid(cellids_end),
            typeid(thrust::counting_iterator<unsigned int>));

  std::vector<double> tmp(4, 0.0);
  Wonton::vector<double> cellval(tmp);
  ASSERT_EQ(typeid(cellval), typeid(thrust::device_vector<double>));

  Wonton::transform(cellids_begin, cellids_end, cellval.begin(),
                    [](int c){return 4.0*c;});
  for (int i = 0; i < 4; i++)
    ASSERT_EQ(4.0*i, cellval[i]);
}

TEST(THRUST_TEST, FOR_EACH) {
  struct aFunctor {
    aFunctor(int size) {vals_.resize(size, 0.0);}
    double operator()(int c) {vals_[c] = 6.0*c;}
    Wonton::vector<double> vals_;
  };

  aFunctor myFunctor(4);

  Wonton::counting_iterator cellids_begin(0);
  Wonton::counting_iterator cellids_end(4);
  Wonton::for_each(cellids_begin, cellids_end, myFunctor);
  for (int i = 0; i < 4; i++)
    ASSERT_EQ(6.0*i, myFunctor.vals_[i]);
}
