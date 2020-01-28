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
  EXPORT wonton-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh/flat
  INCLUDES DESTINATION include/wonton/mesh/flat
  )
      

# HOW DO WE ENSURE THAT JALI IS PROCESSED BEFORE WE ADD THIS TEST

if (ENABLE_UNIT_TESTS)

  if (Jali_DIR)

    # Flat mesh wrapper (has to be paired with a real MPI capable mesh wrapper)
    
    cinch_add_unit(test_flat_mesh_wrapper
      SOURCES flat/test/test_flat_mesh_wrapper.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      # this is relative to the wrapper directory
      INPUTS ../../test_data/shotshell.exo
      ../../test_data/shotshell-v.exo
      POLICY MPI
      THREADS 1)
    
    target_include_directories(test_flat_mesh_wrapper PRIVATE ${Jali_INCLUDE_DIRS})

  endif (Jali_DIR)

endif (ENABLE_UNIT_TESTS)  
