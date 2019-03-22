/*
 This file is part of the Ristra wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_INTPOINT_H_
#define WONTON_INTPOINT_H_

namespace Wonton {

/*!
  @brief An N-dimensional point (of integers rather than doubles)

  A convenient shorthand for collecting indices in N-space.
*/
template<int N>
using IntPoint = std::array<int,N>;

}  // namespace Wonton

#endif  // WONTON_INTPOINT_H_
