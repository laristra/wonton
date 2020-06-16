/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_EQUIFACTOR_H_
#define WONTON_EQUIFACTOR_H_

#include <vector>
#include <algorithm>
#include <unordered_set>
#include <iostream>
#include <cmath>
#include <ctime>
#include <string>

#include "wonton/support/prime_factors.h"

/*!
  @brief Factorize a number N into D equal (or nearly equal) factors
  @param N Number to be factorized
  @param D Number of factors requested
  @param randseed Optional seed for randomization (for reproducibility)
  @returns vector of factors

  The routine is written generally as a bin packing algorithm, but it
  is meant to be used for a regular partitioning of a structured
  mesh/grid. Typically N is the number of processors, and D is the
  number of dimensions along which the mesh is to be partitioned. If N
  can only be decomposed into fewer than D prime factors, some of the
  factors will be 1. The algorithm is like simulated annealing (swap
  random numbers between random min and max sets, allow objective
  function to rise if stalled)

  For the designed usage (parallel partitioning), t is inconceivable
  that we will generate factors greater than what is accommodated by
  an int. If it is, we have to make the 'products' vector and the
  return vector int64_t
*/

namespace Wonton {

#ifdef ENABLE_DEBUG
void print_sets(std::vector<std::vector<int>> sets);
#endif

// gcc 7.3.0 doesn't recognize that this function is used in the code 
// so use the compiler specific attribute to turn off the warning (since we
// use -Werror and cannot get past compilation)
#if defined(__GNUC__) && (__GNUC__ == 7 && __GNUC_MINOR__ == 3)
__attribute__ ((unused))
#endif
static inline
std::vector<int> equifactor(int const N, int const D, int const randseed = 0) {
  clock_t startclock, curclock;
  startclock = clock();

  assert(N > 0 && D > 0);

  if (D == 1)
    return std::vector<int>(1, N);

  // Find all the prime factors of N (with duplicates)

  std::vector<int> prime_facs = prime_factors(N);
  std::sort(prime_facs.begin(), prime_facs.end());
  int const prime_facs_size = prime_facs.size();
  if (prime_facs_size < D) {
    for (int i = prime_facs_size; i < D; i++)
      prime_facs.push_back(1);
  }


  std::vector<std::vector<int>> sets(D, {1});  // init to 1 for degenerate cases
  std::vector<std::vector<int>> minsets;

  // Partition these factors into D sets in a simplistic manner.
  // Since prime_facs is sorted, we are distributing small and large
  // factors equally amongst the sets as much as possible

  int kset = 0;
  for (int const n : prime_facs) {
    sets[kset%D].push_back(n);
    kset++;
  }

#ifdef ENABLE_DEBUG
  // print out initial sets
  std::cerr << "\n\nNumber of sets D " << D << "\n";
  std::cerr << "INITIAL STATE:\n";
  print_sets(sets);
#endif


  std::vector<int> products(D);

  bool outer_done = false;
  int maxiter = 100*D;
  int outer_iter = 0;
  int num_no_change = 0;
  bool allow_swap_climb = false, allow_move_climb = false;
  int nclimbs = 0;
  int olddiff = 0;

  minsets = sets;

  if (randseed)
    srand(randseed);

  while (!outer_done) {

    // compute the measure (product of elements of set) for each set
    products.assign(D, 1.0);
    for (int i = 0; i < D; i++)
      for (int const n : sets[i])
        products[i] *= n;

    // Min and max products
    int minprod = *(std::min_element(products.begin(), products.end()));
    int maxprod = *(std::max_element(products.begin(), products.end()));

    int diff = maxprod-minprod;
    if (diff == 0) {  // all sets have same measure
      outer_done = true;
      continue;
    }

    // collect _all_ sets with min measure and _all_ sets with max measure
    std::vector<int> iminsets, imaxsets;
    for (int i = 0; i < D; i++) {
      if (products[i] == minprod) iminsets.push_back(i);
      if (products[i] == maxprod) imaxsets.push_back(i);
    }

    // Find one random set from the minsets and one from the maxsets
    if (!randseed) {
      curclock = clock();
      srand(curclock - startclock);
    }

    int iset1 = iminsets[rand() % iminsets.size()];
    int iset2 = imaxsets[rand() % imaxsets.size()];

    std::vector<int>& set1 = sets[iset1];
    std::vector<int>& set2 = sets[iset2];
    int & prod1 = products[iset1];
    int & prod2 = products[iset2];


    // swap/move elements between sets to try to reduce difference
    int inner_iter = 0;
    bool inner_done = false;
    while (!inner_done) {

      if (!randseed) {
        curclock = clock();
        srand(curclock - startclock);
      }

      int offset = rand() % set1.size();
      auto j_set1 = set1.begin() + offset;
      if (prod1 != 1) {
        while (*j_set1 == 1) {  // Search for value that's not 1
          offset = rand() % set1.size();
          j_set1 = set1.begin() + offset;
        }
      }
      int& val_set1 = *j_set1;

      if (!randseed) {
        curclock = clock();
        srand(curclock - startclock);
      }

      offset = rand() % set2.size();
      auto j_set2 = set2.begin() + offset;
      if (prod2 != 1) {
        while (*j_set2 == 1) {  // Search for value that's not 1
          offset = rand() % set2.size();
          j_set2 = set2.begin() + offset;
        }
      }
      int& val_set2 = *j_set2;

      // will swapping the two result in a smaller difference in the
      // set measures (product of the two members) or are we giving a
      // random kick because we are stalled?


      int prod1_new = val_set2*prod1/val_set1;
      int prod2_new = val_set1*prod2/val_set2;

      int diff1 = prod2_new - prod1_new;

      if (std::abs(diff1) < std::abs(diff) || allow_swap_climb) {

        std::swap(val_set1, val_set2);
        prod1 = prod1_new;
        prod2 = prod2_new;
        diff = diff1;

        allow_swap_climb = false;

      } else {

        // will moving an element of the max set to the min set reduce
        // the difference or are we giving a random kick because we
        // are stalled?

        prod1_new = prod1*val_set2;
        prod2_new = prod2/val_set2;

        diff1 = prod2_new - prod1_new;

        if (std::abs(diff1) < std::abs(diff) || allow_move_climb) {

          set1.push_back(val_set2);
          val_set2 = 1.0;  // erase is expensive; do a "virtual deletion"
          prod1 = prod1_new;
          prod2 = prod2_new;
          diff = diff1;

          allow_move_climb = false;
        }
      }

      inner_iter++;
      if (unsigned(inner_iter) > 2*set1.size()) inner_done = true;
    }  // while (!inner_done)

    outer_iter++;


    minprod = *(std::min_element(products.begin(), products.end()));
    maxprod = *(std::max_element(products.begin(), products.end()));

    int maxdiff = maxprod-minprod;

    if (maxdiff == olddiff)
      num_no_change++;
    else {
      num_no_change = 0;
      if (maxdiff < olddiff)
        minsets = sets;      // record sets that minimized the min/max diff
      olddiff = maxdiff;
    }

    if (num_no_change > 5) {  // stalled; allow obj. func. to increase
      if (num_no_change%2)
        allow_swap_climb = true;
      else
        allow_move_climb = true;
      nclimbs++;
    }

    if (nclimbs > 10) outer_done = true;
    if (outer_iter > maxiter) outer_done = true;
  }  // while (!outer_done)

#ifdef ENABLE_DEBUG
  std::cerr << "\nFINAL STATE:\n";
  std::vector<std::vector<int>> newsets(D);
  for (int i = 0; i < D; i++) {
    if (minsets[i].size() == 1)
      newsets[i].push_back(minsets[i][0]);
    else
      for (int const n : minsets[i])
        if (n != 1) newsets[i].push_back(n);
  }
  std::swap(minsets,newsets);
  print_sets(minsets);
#endif

  products.assign(D, 1.0);
  for (int i = 0; i < D; i++) {
    products[i] = 1.0;
    for (int const n : minsets[i])
      products[i] *= n;
  }

  return products;
}  // equipartition


#ifdef ENABLE_DEBUG
void print_sets(std::vector<std::vector<int>> sets) {
  int nsets = sets.size();
  int minprod = std::numeric_limits<int>::max();
  int maxprod = std::numeric_limits<int>::min();
  for (int i = 0; i < nsets; i++) {
    std::cerr << "Set " << i << ": ";
    int prod = 1;
    for (int const n : sets[i]) {
      prod *= n;
      std::cerr << n << " ";
    }
    std::cerr << "prod=" << prod << "\n";
    if (prod < minprod)
      minprod = prod;
    if (prod > maxprod)
      maxprod = prod;
  }
  int maxdiff = maxprod - minprod;
  std::cerr << "Max diff: " << maxdiff << "\n";
}
#endif

}  // namespace Wonton

#endif
