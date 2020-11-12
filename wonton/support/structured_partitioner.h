/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_STRUCTURED_PARTITIONER_H_
#define WONTON_STRUCTURED_PARTITIONER_H_


#include "wonton/support/equifactor.h"

namespace Wonton {


/// @brief Partition a D-dimensional structured grid along d axes into N subdomains
/// @tparam D      Mesh/Grid dimensionality
/// @param  N      Number of partitions
/// @param  ncells Number of cells in each axis direction
/// @param  d      Number of directions in which to partition (d <= D)
/// @param  seed   Seed for randomization used in equifactor

template<int D>
std::vector<std::array<std::array<int64_t,D>,2>>
structured_partitioner(int N, std::array<int64_t, D> ncells, int d=D, int seed=0) {

  // factor the number of processors/partitions into D, nearly equal numbers
  std::vector<int> nbins = equifactor(N, d, seed);
  for (int i = d; i < D; i++)
    nbins.push_back(1);

  // compute num cells in each bin along each axis
  std::array<std::vector<int64_t>, D> bincounts;

  for (int i = 0; i < D; i++) {
    int nperbin = ncells[i]/nbins[i];
    int remainder = ncells[i]%nbins[i];

    // spread the cells as evenly as possible
    bincounts[i].assign(nbins[i], nperbin);

    // distribute the remainder as best as possible
    for (int j = 0; j < remainder; j++)
      bincounts[i][j]++;
  }

  // cell offsets of each bin along each axis
  std::array<std::vector<int64_t>, D> binoffsets;

  for (int i = 0; i < D; i++) {
    binoffsets[i].resize(nbins[i]+1);
    binoffsets[i][0] = 0;

    for (int j = 1; j < nbins[i]+1; j++)
      binoffsets[i][j] = binoffsets[i][j-1] + bincounts[i][j-1];
  }

  std::vector<std::array<std::array<int64_t,D>,2>> ilimits(N);

  // One could make this a recursive function but it becomes less readable
  if (D == 1) {
    for (int n = 0; n < N; n++) {
      ilimits[n][0][0] = binoffsets[0][n];
      ilimits[n][1][0] = binoffsets[0][n+1];
    }
  } else if (D == 2) {
    for (int i = 0; i < nbins[0]; i++) {
      for (int j = 0; j < nbins[1]; j++) {
        int n = i*nbins[1]+j;
        ilimits[n][0][0] = binoffsets[0][i];
        ilimits[n][0][1] = binoffsets[1][j];
        ilimits[n][1][0] = binoffsets[0][i+1];
        ilimits[n][1][1] = binoffsets[1][j+1];
      }
    }
  } else if (D == 3) {
    for (int i = 0; i < nbins[0]; i++) {
      int n1 = i*nbins[1];
      for (int j = 0; j < nbins[1]; j++) {
        int n2 = (n1 + j)*nbins[2];
        for (int k = 0; k < nbins[2]; k++) {
          int n = n2 + k;        
          ilimits[n][0][0] = binoffsets[0][i];
          ilimits[n][0][1] = binoffsets[1][j];
          ilimits[n][0][2] = binoffsets[2][k];
          ilimits[n][1][0] = binoffsets[0][i+1];
          ilimits[n][1][1] = binoffsets[1][j+1];
          ilimits[n][1][2] = binoffsets[2][k+1];
        }
      }
    }
  }

  return ilimits;
}

}  // namespace Wonton

#endif
