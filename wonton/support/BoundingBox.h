/*
 This file is part of the Ristra wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_BOUNDING_BOX_H_
#define WONTON_BOUNDING_BOX_H_

<<<<<<< HEAD
=======
#include <array>

>>>>>>> v1.1.2
namespace Wonton {

/*!
  @brief Bounding Box type

  Represents the bounding box of a cell by listing the lower and upper corners
  of the box.
*/
template<int D>
using BoundingBox = std::array<std::array<double,2>,D>;

const int LO = 0;
const int HI = 1;

}  // namespace Wonton

#endif  // WONTON_BOUNDING_BOX_H_
