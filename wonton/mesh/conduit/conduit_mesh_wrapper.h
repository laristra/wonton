/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

#ifndef CONDUIT_MESH_WRAPPER_H_
#define CONDUIT_MESH_WRAPPER_H_

#include <cassert>
#include <algorithm>
#include <vector>
#include <array>
#include <utility>

#include "conduit/conduit.hpp"                  // Conduit header

#include "wonton/mesh/AuxMeshTopology.h"
#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"

namespace Wonton {

std::string conduit_about() {
    return conduit::about();
}

}  // end namespace Wonton

#endif // CONDUIT_MESH_WRAPPER_H_
