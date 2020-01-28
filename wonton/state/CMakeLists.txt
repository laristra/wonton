#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


set(state_HEADERS
  state_vector_base.h
  state_vector_uni.h
  state_vector_multi.h
  state_vector_uni_raw.h
  state_vector_multi_raw.h
  state_manager.h
  )

add_library(wonton_state INTERFACE)

target_include_directories(wonton_state INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include/wonton/state>)


if (Jali_DIR)
  add_subdirectory(jali)
endif (Jali_DIR)

if (ENABLE_FleCSI)
  add_subdirectory(flecsi)
endif ()

add_subdirectory(simple)

add_subdirectory(flat)


install(TARGETS wonton_state
  EXPORT wonton-targets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/state
  INCLUDES DESTINATION include/wonton/state
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_state_vector_uni
    SOURCES test/test_state_vector_uni.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_multi
    SOURCES test/test_state_vector_multi.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_manager
    SOURCES test/test_state_manager.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_uni_raw
    SOURCES test/test_state_vector_uni_raw.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_vector_multi_raw
    SOURCES test/test_state_vector_multi_raw.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_state_manager_raw
    SOURCES test/test_state_manager_raw.cc
    POLICY SERIAL
    )

  # Simple State
  # cinch_add_unit(test_simple_state
  #   SOURCES simple/test/test_simple_state.cc
  #   POLICY SERIAL
  #   )
  
  # # Simple State wrapper
  # cinch_add_unit(test_simple_state_wrapper
  #   SOURCES simple/test/test_simple_state_wrapper.cc
  #   POLICY SERIAL)
  
  # # Simple State MMwrapper
  # cinch_add_unit(test_simple_state_mm_wrapper
  #   SOURCES simple/test/test_simple_state_mm_wrapper.cc
  #   POLICY SERIAL)
  
  
#  if (Jali_DIR)
#    cinch_add_unit(test_jali_state_wrapper
#      SOURCES jali/test/test_jali_state_wrapper.cc
#      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
#      POLICY MPI
#      THREADS 1)
    
    # Flat state wrapper
#    cinch_add_unit(test_flat_state_wrapper
#      SOURCES flat/test/test_flat_state_wrapper.cc
#      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
#      POLICY MPI
#      THREADS 1)
#  endif (Jali_DIR)
  
  
  # HAVE TO ADD FLECSI STATE WRAPPER TESTS
  
endif (ENABLE_UNIT_TESTS)
