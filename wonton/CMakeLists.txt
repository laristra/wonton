#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_subdirectory(support)

if (WONTON_ENABLE_MPI)
  add_subdirectory(distributed)
endif ()

add_subdirectory(intersect)
add_subdirectory(mesh)
add_subdirectory(state)
add_subdirectory(swarm)
