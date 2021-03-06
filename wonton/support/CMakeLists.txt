#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(support)

add_library(wonton_support)

target_sources(wonton_support PRIVATE Matrix.cc svd.cc)

target_include_directories(wonton_support PUBLIC
  $<BUILD_INTERFACE:${wonton_SOURCE_DIR}>
  $<BUILD_INTERFACE:${wonton_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(wonton_support_HEADERS
  BoundingBox.h
  CoordinateSystem.h
  equifactor.h
  lsfits.h
  Matrix.h
  moment_index.h
  Point.h
  Polytope.h
  prime_factors.h
  structured_partitioner.h
  svd.h
  Vector.h
  wonton.h
)
set_target_properties(wonton_support PROPERTIES
  PUBLIC_HEADER "${wonton_support_HEADERS}")

if (WONTON_ENABLE_MPI)
  target_link_libraries(wonton_support PUBLIC MPI::MPI_CXX)

  # Need this compile definition for pulling in correct parts of the code
  target_compile_definitions(wonton_support PUBLIC WONTON_ENABLE_MPI)
  
  # MPI flag to avoid pulling in MPI_CXX libraries (Moot after 3.0.0)
  target_compile_definitions(wonton_support PUBLIC OMPI_SKIP_MPICXX)
endif ()

#-----------------------------------------------------------------------------
# Lapacke
#-----------------------------------------------------------------------------
if (WONTON_ENABLE_LAPACKE AND LAPACKE_FOUND)
  target_include_directories(wonton_support PUBLIC ${LAPACKE_INCLUDE_DIRS})
  target_compile_definitions(wonton_support PUBLIC WONTON_HAS_LAPACKE)
  target_link_libraries(wonton_support PUBLIC ${LAPACKE_LIBRARIES})
endif ()

#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
if (WONTON_ENABLE_THRUST)   # if it is overridden by the command line

  target_link_libraries(wonton_support PUBLIC ${THRUST_LIBRARIES})
  target_include_directories(wonton_support PUBLIC ${THRUST_INCLUDE_DIRS})

  # Need this compile definition for pulling in correct parts of the code
  target_compile_definitions(wonton_support PUBLIC WONTON_ENABLE_THRUST)

  target_compile_definitions(wonton_support PUBLIC THRUST_HOST_SYSTEM=${THRUST_HOST_BACKEND})
  target_compile_definitions(wonton_support PUBLIC THRUST_DEVICE_SYSTEM=${THRUST_DEVICE_BACKEND})

else ()

  if (Boost_LIBRARIES)
    target_link_libraries(wonton_support PUBLIC ${Boost_LIBRARIES})
  else ()
    target_include_directories(wonton_support PUBLIC ${Boost_INCLUDE_DIRS})
  endif ()

endif()

# TCMalloc is useful in general but particularly with multi-core runs
if (WONTON_ENABLE_TCMALLOC)
  target_link_libraries(wonton_support PUBLIC ${TCMALLOC_LIBRARIES})
endif ()

#----------------------------------------------------------------------------
# Kokkos dependencies
#----------------------------------------------------------------------------

if (WONTON_ENABLE_Kokkos)
  target_link_libraries(wonton_support PUBLIC ${Kokkos_LIBRARIES})

  # Do we need this? Shouldn't Kokkos take care of this?
  if (WONTON_ENABLE_CUDA)
    target_compile_definitions(wonton_support PUBLIC ${CUDA_NVCC_FLAGS})
  endif ()
endif ()


target_link_libraries(wonton INTERFACE wonton_support)

install(TARGETS wonton_support
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/support
  INCLUDES DESTINATION include/wonton/support
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  wonton_add_unittest(test_point
    SOURCES test/test_point.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  wonton_add_unittest(test_vector
    SOURCES test/test_vector.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  wonton_add_unittest(test_matfuncs
    SOURCES test/test_matfuncs.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  wonton_add_unittest(test_coord_sys
    SOURCES test/test_coord_sys.cc
    LIBRARIES wonton_support
    POLICY SERIAL
    )

  wonton_add_unittest(test_moment_index
    SOURCES test/test_moment_index.cc
    LIBRARIES wonton_support
    POLICY SERIAL
    )

  wonton_add_unittest(test_polytope_2
    SOURCES test/test_polytope_2.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  wonton_add_unittest(test_polytope_3
    SOURCES test/test_polytope_3.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  if (WONTON_ENABLE_THRUST)
    wonton_add_unittest(test_thrust
      SOURCES test/test_thrust.cc
      LIBRARIES wonton_support
      POLICY SERIAL)
  endif ()

  wonton_add_unittest(test_prime_factors
    SOURCES test/test_prime_factors.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  wonton_add_unittest(test_equifactor
    SOURCES test/test_equifactor.cc
    LIBRARIES wonton_support
    POLICY SERIAL)
  
  wonton_add_unittest(test_structured_partitioner
    SOURCES test/test_structured_partitioner.cc
    LIBRARIES wonton_support
    POLICY SERIAL)

  if (WONTON_ENABLE_MPI)
    wonton_add_unittest(test_structured_partitioner_par
      SOURCES test/test_structured_partitioner_par.cc
      LIBRARIES wonton_support
      POLICY MPI
      THREADS 4)
  endif ()
  
endif(ENABLE_UNIT_TESTS)
