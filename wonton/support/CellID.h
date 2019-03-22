/*
 This file is part of the Ristra wonton project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef WONTON_CELLID_H_
#define WONTON_CELLID_H_

namespace Wonton {

/*!
  @brief Cell ID type

  Makes it easier if we change the type used for cell IDs in the future.  Also
  helps clarify that this particular int64_t isn't just a random integer, but
  is actually a cell ID.
*/
using CellID = int64_t;

}  // namespace Wonton

#endif  // WONTON_CELLID_H_
