#[[
This file is part of the Ristra Wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]
#-----------------------------------------------------------------------------~#
project(mesh)


# Only header files in this directory and it's subdirectories
set(mesh_HEADERS
  AuxMeshTopology.h
  )


add_library(wonton_mesh INTERFACE)

target_include_directories(wonton_mesh INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)


if (Jali_DIR)
  add_subdirectory(jali)
endif (Jali_DIR)

if (FleCSI_DIR)
  add_subdirectory(flecsi)
endif (FleCSI_DIR)

add_subdirectory(simple)

add_subdirectory(adaptive_refinement)

add_subdirectory(direct_product)

target_link_libraries(wonton INTERFACE wonton_mesh)

install(TARGETS wonton_mesh
  EXPORT wonton_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/wonton/mesh
  INCLUDES DESTINATION include/wonton/mesh
  )

if (ENABLE_UNIT_TESTS)

  # Adaptive Refinement Mesh
  # cinch_add_unit(test_adaptive_refinement_mesh
  #   SOURCES adaptive_refinement/test/test_adaptive_refinement_mesh.cc
  #   LIBRARIES wonton_mesh
  #   POLICY SERIAL
  #   )

  # Adaptive Refinement Mesh Wrapper
  # cinch_add_unit(test_adaptive_refinement_mesh_wrapper
  #   SOURCES adaptive_refinement/test/test_adaptive_refinement_mesh_wrapper.cc
  #   LIBRARIES wonton_mesh
  #   POLICY SERIAL
  #   )

  # Direct Product Mesh
#   cinch_add_unit(test_direct_product_mesh
#     SOURCES direct_product/test/test_direct_product_mesh.cc
#     LIBRARIES wonton_mesh
#     POLICY SERIAL
#     )

# # Direct Product Mesh Wrapper
#   cinch_add_unit(test_direct_product_mesh_wrapper
#     SOURCES direct_product/test/test_direct_product_mesh_wrapper.cc
#     LIBRARIES wonton_mesh
#     POLICY SERIAL
#     )

  # Simple Mesh
  # cinch_add_unit(test_simple_mesh
  #   SOURCES simple/test/test_simple_mesh.cc
  #   LIBRARIES wonton_mesh
  #   POLICY SERIAL
  #   )

  # Simple Mesh wrapper
  # cinch_add_unit(test_simple_mesh_wrapper
  #   SOURCES simple/test/test_simple_mesh_wrapper.cc
  #   LIBRARIES wonton_mesh
  #   POLICY SERIAL)

  if (Jali_DIR)
    # cinch_add_unit(test_jali_mesh_wrapper
    #   SOURCES jali/test/test_jali_mesh_wrapper.cc
    #   LIBRARIES wonton_mesh
    #   # this is relative to the wrapper directory
    #   INPUTS ../../test_data/shotshell.exo
    #   ../../test_data/shotshell-v.exo
    #   POLICY MPI
    #   THREADS 1)

    # Flat mesh wrapper (has to be paired with a real MPI capable mesh
    # wrapper)
    # cinch_add_unit(test_flat_mesh_wrapper
    #   SOURCES flat/test/test_flat_mesh_wrapper.cc
    #   LIBRARIES wonton_mesh
    #   # this is relative to the wrapper directory
    #   INPUTS ../../test_data/shotshell.exo
    #   ../../test_data/shotshell-v.exo
    #   POLICY MPI
    #   THREADS 1)
  endif (Jali_DIR)

  if (ENABLE_FleCSI)
    cinch_add_unit(test_flecsi_mesh_wrapper
      SOURCES flecsi/test/test_flecsi_mesh_wrapper.cc
      LIBRARIES wonton_mesh
      POLICY MPI
      THREADS 1)
  endif (ENABLE_FleCSI)
endif(ENABLE_UNIT_TESTS)
