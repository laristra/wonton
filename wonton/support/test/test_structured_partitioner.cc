/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include <iostream>

#include "gtest/gtest.h"
#include "wonton/support/structured_partitioner.h"

TEST(StructPartitioner, Serial) {

  {
    // Partition 2x2 mesh into 4 partitions (trivial)
    std::array<int, 2> ncells = {2, 2};
    auto partlimits = Wonton::structured_partitioner<4, 2>(ncells);
    
    std::array<std::array<int, 4>, 4> exp_partlimits =
        { {{0,0,1,1}, {0,1,1,2}, {1,0,2,1}, {1,1,2,2}} };
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
           ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }

  {
    // Partition 2x3 mesh into 4 partitions
    std::array<int, 2> ncells = {2, 3};

    auto partlimits = Wonton::structured_partitioner<4, 2>(ncells);

    std::array<std::array<int, 4>, 4> exp_partlimits =
        { {{0,0,1,2}, {0,2,1,3}, {1,0,2,2}, {1,2,2,3}} };
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }

  {
    // Partition 55x10 mesh into 3 partitions
    std::array<int, 2> ncells = {55, 10};

    auto partlimits = Wonton::structured_partitioner<3, 2>(ncells);

    std::array<std::array<int, 4>, 3> exp_partlimits =
        { {{0,0,19,10}, {19,0,37,10}, {37,0,55,10}} };  // Removing the outermost braces causes compilation error
    
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 4; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }
  
  {
    // Partition 55x10 mesh into 4 partitions BUT ONLY IN 1 DIR
    std::array<int, 2> ncells = {55, 10};

    auto partlimits = Wonton::structured_partitioner<4, 2, 1>(ncells);

    std::array<std::array<int, 4>, 4> exp_partlimits =
        { {{0,0,14,10}, {14,0,28,10}, {28,0,42,10}, {42,0,55,10}} };
    
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }


  {
    // Partition 55x10x12 mesh into 8 partitions
    std::array<int, 3> ncells = {55, 10, 12};

    auto partlimits = Wonton::structured_partitioner<8, 3>(ncells);

    std::array<std::array<int, 6>, 8> exp_partlimits =
        { {{0,0,0,28,5,6}, {0,0,6,28,5,12},
           {0,5,0,28,10,6}, {0,5,6,28,10,12},
           {28,0,0,55,5,6}, {28,0,6,55,5,12},
           {28,5,0,55,10,6}, {28,5,6,55,10,12}} };
           
    
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 6; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }


  {
    // Partition 55x10x12 mesh into 8 partitions BUT ONLY IN 2 DIRS
    std::array<int, 3> ncells = {55, 10, 12};

    auto partlimits = Wonton::structured_partitioner<8, 3, 2>(ncells);

    std::array<std::array<int, 6>, 8> exp_partlimits =
        { {{0,0,0,14,5,12}, {0,5,0,14,10,12},
           {14,0,0,28,5,12}, {14,5,0,28,10,12},
           {28,0,0,42,5,12}, {28,5,0,42,10,12},
           {42,0,0,55,5,12}, {42,5,0,55,10,12}} };
           
    
    for (int i = 0; i < 8; i++)
      for (int j = 0; j < 6; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }


  {
    // Partition 55x10x5 mesh into 5 partitions BUT ONLY IN 1 DIR
    std::array<int, 3> ncells = {55, 10, 5};

    auto partlimits = Wonton::structured_partitioner<5, 3, 1>(ncells);

    std::array<std::array<int, 6>, 5> exp_partlimits =
        { {{0,0,0,11,10,5}, {11,0,0,22,10,5}, {22,0,0,33,10,5},
           {33,0,0,44,10,5}, {44,0,0,55,10,5}} };
    
    for (int i = 0; i < 5; i++)
      for (int j = 0; j < 6; j++)
        ASSERT_EQ(exp_partlimits[i][j], partlimits[i][j]);
  }
}

