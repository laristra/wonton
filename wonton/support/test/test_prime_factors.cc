/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include <iostream>

#include "gtest/gtest.h"
#include "wonton/support/prime_factors.h"

/*!
  Test if we can get the prime factors (with repetitions) of a number. Tests up to the first 10 prime factors
*/

TEST(PrimeFactors, PrimeFactors) {

  std::vector<int> factors, expected_factors;

  factors = Wonton::prime_factors(0);
  ASSERT_EQ(0, factors.size());

  factors = Wonton::prime_factors(1);
  ASSERT_EQ(1, factors.size());

  factors = Wonton::prime_factors(2);
  ASSERT_EQ(1, factors.size());
  ASSERT_EQ(2, factors[0]);

  factors = Wonton::prime_factors(3);
  ASSERT_EQ(1, factors.size());
  ASSERT_EQ(3, factors[0]);

  factors = Wonton::prime_factors(4);
  expected_factors = std::vector<int>({2,2});
  ASSERT_EQ(expected_factors.size(), factors.size());
  for (int i = 0; i < expected_factors.size(); i++)
    ASSERT_EQ(expected_factors[i], factors[i]);

  factors = Wonton::prime_factors(36);
  expected_factors = std::vector<int>({2,2,3,3}); 
  ASSERT_EQ(expected_factors.size(), factors.size());
  for (int i = 0; i < expected_factors.size(); i++)
    ASSERT_EQ(expected_factors[i], factors[i]);

  factors = Wonton::prime_factors(36);
  expected_factors = std::vector<int>({2,2,3,3}); 
  ASSERT_EQ(expected_factors.size(), factors.size());
  for (int i = 0; i < expected_factors.size(); i++)
    ASSERT_EQ(expected_factors[i], factors[i]);

  factors = Wonton::prime_factors(7350);
  expected_factors = std::vector<int>({2,3,5,5,7,7}); 
  ASSERT_EQ(expected_factors.size(), factors.size());
  for (int i = 0; i < expected_factors.size(); i++)
    ASSERT_EQ(expected_factors[i], factors[i]);

  factors = Wonton::prime_factors(3509);
  expected_factors = std::vector<int>({11,11,29}); 
  ASSERT_EQ(expected_factors.size(), factors.size());
  for (int i = 0; i < expected_factors.size(); i++)
    ASSERT_EQ(expected_factors[i], factors[i]);

}
