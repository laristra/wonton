/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#include "wonton/mesh/conduit/conduit_mesh_wrapper.h"

#include <iostream>

#include "gtest/gtest.h"

#include "wonton/support/Point.h"

/*!
  file test_conduit_mesh_wrapper.cc
  @brief Unit tests for the Conduit mesh wrapper class
 */

TEST(Conduit_Mesh_Wrapper, Sanity_Test) {
  ASSERT_NE(Wonton::conduit_about(), "");
}
