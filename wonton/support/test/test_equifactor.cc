/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include <iostream>

#include "gtest/gtest.h"
#include "wonton/support/equifactor.h"

TEST(EquiFactor, EquiFactor) {
  std::vector<int> factors;
  
  // Try factoring a square
  factors = Wonton::equifactor(64, 2);
  for (int i = 0; i < 2; i++)
    ASSERT_EQ(factors[i], 8);
  
  
  // Try factoring a cube
  factors = Wonton::equifactor(216, 3);
  for (int i = 0; i < 3; i++)
    ASSERT_EQ(factors[i], 6);
  
  // Try a general number with a general number of partitions
  factors = Wonton::equifactor(144, 4);  // should give us 4, 4, 3, 3
  std::sort(factors.begin(), factors.end());
  ASSERT_EQ(factors[0], 3);
  ASSERT_EQ(factors[1], 3);
  ASSERT_EQ(factors[2], 4);
  ASSERT_EQ(factors[3], 4);
  
  // Another general number 10*10*11*11*12*12*13*13 = 294465600
  // Specify the randomization seed to get reproducibility
  int P = 10*10*11*11*12*12*13*13;
  factors = Wonton::equifactor(P, 5, 42);  // (44, 45, 52, 52, 55)
  std::sort(factors.begin(), factors.end());
  ASSERT_EQ(factors[0], 44);
  ASSERT_EQ(factors[1], 45);
  ASSERT_EQ(factors[2], 52);
  ASSERT_EQ(factors[3], 52);
  ASSERT_EQ(factors[4], 55);
  
  // How about a prime number?
  factors = Wonton::equifactor(13, 3);  // 13, 1, 1
  std::sort(factors.begin(), factors.end());
  ASSERT_EQ(factors[0], 1);
  ASSERT_EQ(factors[1], 1);
  ASSERT_EQ(factors[2], 13);
  
}
