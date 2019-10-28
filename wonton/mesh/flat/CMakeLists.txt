#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

set(flat_mesh_wapper_HEADERS
  flat_mesh_wrapper.h
  )

target_sources(mesh INTERFACE ${flat_mesh_wrapper_HEADERS})

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
