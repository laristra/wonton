/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STRUCTURED_PARTITIONER_H
#define WONTON_STRUCTURED_PARTITIONER_H_


// namespace {

// /* given partition offsets (see structured_partitioner below for
//  * meaning of partition offsets) in each direction of a structured
//  * grid, fill in lower and upper limits of each partition limits */

// template<int N, int D>
// void fill_limits(std::array<std::vector<int>, D> const & binoffsets) {}


// template<int N>
// std::array<std::array<int, 2>, N>
// fill_limits<1>(std::array<std::vector<int>, 1> const & binoffsets) {
//   std::array<std::array<int, 2>, N> ilimits;
//   for (int n = 0; n < N; n++) {
//     ilimits[n][0] = binoffsets[n];
//     ilimits[n][1] = binoffsets[n+1];
//   }
//   return ilimits;
// }

// template<int N>
// std::array<std::array<int, 4>, N>
// fill_limits<2>(std::array<std::vector<int>, 2> const & binoffsets) {
//   std::array<int, 2> nbins = {binoffsets[0].size()-1, binoffsets[1].size()-1};
//   std::array<std::array<int, 4>, N> ilimits;

//   for (int i = 0; i < nbins[0]; i++) {
//     for (int j = 0; j < nbins[1]; j++) {
//       int n = i*nbins[1]+j;
//       ilimits[n][0] = binoffsets[0][i];
//       ilimits[n][1] = binoffsets[1][j];
//       ilimits[n][2] = binoffsets[0][i+1];
//       ilimits[n][3] = binoffsets[1][j+1];
//     }
//   }
// }

// template<int N>
// std::array<std::array<int, 6>, N>
// fill_limits<3>(std::array<std::vector<int>, 3> const & binoffsets) {
//   std::array<int, 3> nbins = {binoffsets[0].size()-1, binoffsets[1].size()-1,
//                               binoffsets[2].size()-1};
//   std::array<std::array<int, 6>, N> ilimits;

//   for (int i = 0; i < nbins[0]; i++) {
//     int n1 = i*nbins[1];
//     for (int j = 0; j < nbins[1]; j++) {
//       int n2 = (n1 + j)*nbins[2];
//       for (int k = 0; k < nbins[2]; k++) {
//         int n = n2 + k;        
//         ilimits[n][0] = binoffsets[0][i];
//         ilimits[n][1] = binoffsets[1][j];
//         ilimits[n][2] = binoffsets[2][k]
//         ilimits[n][3] = binoffsets[0][i+1];
//         ilimits[n][4] = binoffsets[1][j+1];
//         ilimits[n][5] = binoffsets[2][k+1];
//       }
//     }
//   }
// }

// }



#include "wonton/support/equifactor.h"

namespace Wonton {


/* Partition a D-dimensional structured grid along d axes into N subdomains */

template<int N, int D, int d=D>
std::array<std::array<int, 2*D>, N>
structured_partitioner(std::array<int, D> ncells, int randseed = 0) {

  std::vector<int> nbins = equifactor(N, d, randseed);
  for (int i = d; i < D; i++)
    nbins.push_back(1);

  std::array<std::vector<int>, D> bincounts;  // num cells in each bin along each axis

  for (int i = 0; i < D; i++) {
    int nperbin = ncells[i]/nbins[i];
    int remainder = ncells[i]%nbins[i];

    // spread the cells as evenly as possible
    bincounts[i].assign(nbins[i], nperbin);

    // distribute the remainder as best as possible
    for (int j = 0; j < remainder; j++)
      bincounts[i][j]++;
  }

  std::array<std::vector<int>, D> binoffsets; // cell offsets of each bin along each axis

  for (int i = 0; i < D; i++) {
    binoffsets[i].resize(nbins[i]+1);
    binoffsets[i][0] = 0;

    for (int j = 1; j < nbins[i]+1; j++)
      binoffsets[i][j] = binoffsets[i][j-1] + bincounts[i][j-1];
  }

  std::array<std::array<int, 2*D>, N> ilimits;

  // One could make this a recursive function but it becomes less readable
  if (D == 1) {
    for (int n = 0; n < N; n++) {
      ilimits[n][0] = binoffsets[0][n];
      ilimits[n][1] = binoffsets[0][n+1];
    }
  } else if (D == 2) {
    for (int i = 0; i < nbins[0]; i++) {
      for (int j = 0; j < nbins[1]; j++) {
        int n = i*nbins[1]+j;
        ilimits[n][0] = binoffsets[0][i];
        ilimits[n][1] = binoffsets[1][j];
        ilimits[n][2] = binoffsets[0][i+1];
        ilimits[n][3] = binoffsets[1][j+1];
      }
    }
  } else if (D == 3) {
    for (int i = 0; i < nbins[0]; i++) {
      int n1 = i*nbins[1];
      for (int j = 0; j < nbins[1]; j++) {
        int n2 = (n1 + j)*nbins[2];
        for (int k = 0; k < nbins[2]; k++) {
          int n = n2 + k;        
          ilimits[n][0] = binoffsets[0][i];
          ilimits[n][1] = binoffsets[1][j];
          ilimits[n][2] = binoffsets[2][k];
          ilimits[n][3] = binoffsets[0][i+1];
          ilimits[n][4] = binoffsets[1][j+1];
          ilimits[n][5] = binoffsets[2][k+1];
        }
      }
    }
  }

  return ilimits;
}

}  // namespace Wonton

#endif
