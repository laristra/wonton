/*
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
*/

// All the include files we want to check to make sure they won't
// cause multiply defined symbols. If some function does cause multiply
// defined symbol errors, the options are to
//
// 0. Make sure you have include guards so that multiple inclusions of
// the include file in a source file is ok
// 1. Move it to a .cc file
// 2. Make it inline
// 3. Make it static
// 4. Enclose it in its own namespace
//
// Option 3 and 4 will cause each translation unit (compiled source file) to
// have it's own copy of the function


#ifndef WONTON_MULTIDEF_H_
#define WONTON_MULTIDEF_H_

#include "wonton/support/wonton.h"
#include "wonton/support/Point.h"
#include "wonton/support/Vector.h"
#include "wonton/support/BoundingBox.h"
#include "wonton/support/CoordinateSystem.h"
#include "wonton/support/Matrix.h"

#include "wonton/mesh/AuxMeshTopology.h"

#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh.h"
#include "wonton/mesh/adaptive_refinement/adaptive_refinement_mesh_wrapper.h"

#include "wonton/mesh/flat/flat_mesh_wrapper.h"

#include "wonton/mesh/direct_product/direct_product_mesh.h"
#include "wonton/mesh/direct_product/direct_product_mesh_wrapper.h"

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

#include "wonton/state/state_vector_base.h"
#include "wonton/state/state_vector_uni.h"
#include "wonton/state/state_vector_multi.h"
#include "wonton/state/state_vector_uni_raw.h"
#include "wonton/state/state_vector_multi_raw.h"
#include "wonton/state/state_manager.h"

#include "wonton/state/flat/flat_state_mm_wrapper.h"

#include "wonton/state/simple/simple_state.h"
#include "wonton/state/simple/simple_state_mm_wrapper.h"

#endif
