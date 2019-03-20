#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (Jali_DIR)
    set(jali_state_wrapper_headers
        jali/jali_state_wrapper.h
        )
endif (Jali_DIR)

if (ENABLE_FleCSI)
    set(flecsi_state_wrapper_headers
        flecsi/flecsi_state_wrapper.h
        )
endif ()

set (state_HEADERS
     state_vector_base.h
     state_vector_uni.h
     state_vector_multi.h
     state_vector_uni_raw.h
     state_vector_multi_raw.h
     state_manager.h
     simple/simple_state.h
     simple/simple_state_wrapper.h
     simple/simple_state_mm_wrapper.h
     flat/flat_state_wrapper.h
     flat/flat_state_mm_wrapper.h
     ${jali_state_wrapper_headers}
     ${flecsi_state_wrapper_headers}
     PARENT_SCOPE
     )

# Add source files
set (wonton_state_SOURCES
  PARENT_SCOPE
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
  cinch_add_unit(test_simple_state
    SOURCES simple/test/test_simple_state.cc
    POLICY SERIAL
    )
  
  # Simple State wrapper
  cinch_add_unit(test_simple_state_wrapper
    SOURCES simple/test/test_simple_state_wrapper.cc
    POLICY SERIAL)
  
  # Simple State MMwrapper
  cinch_add_unit(test_simple_state_mm_wrapper
    SOURCES simple/test/test_simple_state_mm_wrapper.cc
    POLICY SERIAL)
  
  
  if (Jali_DIR)
    cinch_add_unit(test_jali_state_wrapper
      SOURCES jali/test/test_jali_state_wrapper.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
    
    # Flat state wrapper
    cinch_add_unit(test_flat_state_wrapper
      SOURCES flat/test/test_flat_state_wrapper.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1)
  endif (Jali_DIR)
  
  
  # HAVE TO ADD FLECSI STATE WRAPPER TESTS
  
endif (ENABLE_UNIT_TESTS)
