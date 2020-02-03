#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

set(flat_mesh_wapper_HEADERS
  flat_mesh_wrapper.h
  )

add_library(flat_mesh_wrapper INTERFACE)


target_include_directories(flat_mesh_wrapper INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

# add this header-only library as a dependency of upper level mesh library

target_link_libraries(wonton_mesh INTERFACE flat_mesh_wrapper)

install(TARGETS flat_state_wrapper
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/flat
  INCLUDES DESTINATION include/wonton/mesh/flat
  )
      

if (ENABLE_UNIT_TESTS)

  if (ENABLE_Jali)

    # Flat mesh wrapper (has to be paired with a real MPI capable mesh wrapper)
    
    cinch_add_unit(test_flat_mesh_wrapper
      SOURCES test/test_flat_mesh_wrapper.cc
      LIBRARIES flat_mesh_wrapper jali_mesh_wrapper
      # this is relative to the wrapper directory
      INPUTS ${CMAKE_SOURCE_DIR}/test_data/shotshell.exo
      ${CMAKE_SOURCE_DIR}/test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)
    
  endif (ENABLE_Jali)

endif (ENABLE_UNIT_TESTS)  
