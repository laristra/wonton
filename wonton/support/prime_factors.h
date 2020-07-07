/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_PRIME_FACTOR_H_
#define WONTON_PRIME_FACTOR_H_

namespace Wonton {

/*!
  @brief Compute the prime factors of number N (including repeated factors) 
  @param N     Number to be factorized
  @returns     vector of factors (with repeated factors)
*/

static
std::vector<int> prime_factors(unsigned int const N) {
  constexpr int nprimes = 10;
  int primes[nprimes] = {2,3,5,7,11,13,17,19,23,29};  // unlikely to need more
  std::vector<int> factors;
  
  switch (N) {
    case 0: case 1:
      break;
    case 2: case 3:
      factors.push_back(N);
      break;
    default: {
      int num = N;
      int i = 0;
      bool done = false;
      while (!done) {
        if (num%primes[i] == 0) {
          factors.push_back(primes[i]);
          num = num/primes[i];
        } else
          i = i+1;
        if (num == 1 || i == nprimes)
          done = true;
      }
      break;
    }
  }

  return factors;
}
  
}  // end namespace Wonton

#endif
