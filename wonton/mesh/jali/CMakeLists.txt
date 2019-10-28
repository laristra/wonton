#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#

if (Jali_DIR)
  set(jali_mesh_wrapper_HEADERS
    jali_mesh_wrapper.h
    )

  target_sources(mesh INTERFACE ${jali_mesh_wrapper_HEADERS})

  # Jali dependencies
  target_include_directories(mesh INTERFACE ${Jali_INCLUDE_DIRS})
  target_link_libraries(mesh INTERFACE ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})

  cinch_add_unit(test_jali_mesh_wrapper
    SOURCES test/test_jali_mesh_wrapper.cc
    LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
    # this is relative to the wrapper directory
    INPUTS ../../../test_data/shotshell.exo
    ../../../test_data/shotshell-v.exo
    POLICY MPI
    THREADS 1)

  # unit test needs Jali includes to compile
  target_include_directories(test_jali_mesh_wrapper PRIVATE ${Jali_INCLUDE_DIRS})

endif (Jali_DIR)
