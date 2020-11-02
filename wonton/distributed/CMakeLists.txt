#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(wonton_distributed)

add_library(wonton_distributed INTERFACE)

target_include_directories(wonton_distributed INTERFACE
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/wonton/distributed>)

set(wonton_distributed_HEADERS
  mpi_ghost_manager.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(wonton_distributed PROPERTIES
#   PUBLIC_HEADER ${wonton_distributed_HEADERS}
#   )
#
# Directly install files instead
install(FILES ${wonton_distributed_HEADERS} DESTINATION include/wonton/distributed)

target_link_libraries(wonton_distributed INTERFACE wonton_mesh wonton_state)

target_link_libraries(wonton INTERFACE wonton_distributed)

install(TARGETS wonton_distributed
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/distributed
  INCLUDES DESTINATION include/wonton/distributed
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  if (WONTON_ENABLE_Jali)
    wonton_add_unittest(test_mpi_ghost_manager
      SOURCES test/test_mpi_ghost_manager.cc
      LIBRARIES wonton_distributed
      POLICY MPI
      THREADS 4
      )
  endif ()
  
endif ()


