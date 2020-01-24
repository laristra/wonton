/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include <iostream>

#include "gtest/gtest.h"
#include "wonton/support/wonton.h"
#include "wonton/support/structured_partitioner.h"

TEST(StructuredPartioner, Parallel) {
  {
    // parallel test
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (nprocs > 1) {
      // Partition 101x321x529 cell mesh into 16 partitions on every
      // processor and make sure we got the exact same result (the
      // number of processors we run this on has no relation to how
      // many partitions we are requesting).

      int seed = 42;
      std::array<int, 3> ncells = {144, 256, 16};
      auto partlimits = Wonton::structured_partitioner<3>(16, ncells, 3, seed);

      int limitsarray[96];  // flattened array
      for (int i = 0; i < 16; i++)
        for (int j = 0; j < 3; j++) {
          limitsarray[i*6+j]   = partlimits[i][0][j];  // lower limit
          limitsarray[i*6+j+3] = partlimits[i][1][j];  // upper limit
        }

      int limitsarray_all[384];
      MPI_Allgather(limitsarray, 96, MPI_INT, limitsarray_all, 96, MPI_INT,
                    MPI_COMM_WORLD);

      for (int i = 0; i < 16; i++)
        for (int j = 0; j < 6; j++)
          for (int p = 1; p < nprocs; p++)
            ASSERT_EQ(limitsarray_all[i*6+j], limitsarray_all[p*96+i*6+j]);
    }
  }
}
