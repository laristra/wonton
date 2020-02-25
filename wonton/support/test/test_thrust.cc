/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/


#include "gtest/gtest.h"

#include "wonton/support/wonton.h"
#include <iostream>

TEST(THRUST_TEST, TRANSFORM) {
  Wonton::counting_iterator cellids_begin(0);
  Wonton::counting_iterator cellids_end(4);
  ASSERT_EQ(typeid(cellids_begin),
            typeid(thrust::counting_iterator<unsigned int>));
  ASSERT_EQ(typeid(cellids_end),
            typeid(thrust::counting_iterator<unsigned int>));

  Wonton::vector<double> cellval(4, 0.0);
  ASSERT_EQ(typeid(cellval), typeid(thrust::device_vector<double>));

  Wonton::transform(cellids_begin, cellids_end, cellval.begin(),
                    [](int c){return 4.0*c;});
  for (int i = 0; i < 4; i++)
    ASSERT_EQ(4.0*i, cellval[i]);
}

TEST(THRUST_TEST, FOR_EACH) {
  Wonton::vector<double> cellvals(4, 0.0);
  auto multby6 = [&cellvals](int c) {cellvals[c] = 6.0*c;};

  Wonton::counting_iterator cellids_begin(0);
  Wonton::counting_iterator cellids_end(4);

  std::for_each(cellids_begin, cellids_end, multby6);

  for (int i = 0; i < 4; i++)
    ASSERT_EQ(6.0*i, cellvals[i]);
}
