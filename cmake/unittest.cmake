#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


#----------------------------------------------------------------------------
# Google Test
#----------------------------------------------------------------------------

find_package(GTest QUIET)  # This will catch externally installed GTest

if (NOT GTest_FOUND OR NOT TARGET GTest::gtest)  # build from submodule
  find_package(Threads)  # Find thread libraries for system
  
  add_subdirectory(googletest)
endif ()

set(GTEST_LIBRARIES GTest::gtest)

#[===========================================================================[
.. command:: wonton_add_unit

The ``wonton_add_unit`` function creates a custom unit test with
various runtime policies::

wonton_add_unit(<name> [<option>...])

General options are:

``SOURCES <sources>...``
The sources necessary to build the test executable
``INPUTS <inputs>...``
The input files used to run the test
``POLICY <policy>``
The runtime policy to use when executing the test (SERIAL, MPI)
``THREADS <threads>...``
The number of threads to run the test with
``LIBRARIES <libraries>...``
List of libraries to link target against
``DEFINES <defines>...``
Defines to set when building test target
``ARGUMENTS <test arguments>
Arguments supplied to the test command line
#]===========================================================================]

function(wonton_add_unit name)

  #--------------------------------------------------------------------------#
  # Setup argument options.
  #--------------------------------------------------------------------------#
  
  set(one_value_args POLICY)
  set(multi_value_args SOURCES INPUTS LIBRARIES DEFINES THREADS)
  cmake_parse_arguments(unit "${one_value_args}" "${multi_value_args}" ${ARGN})
  
  #--------------------------------------------------------------------------#
  # Make sure that the user specified sources and the list contains Main.cc
  #--------------------------------------------------------------------------#
  
  if (NOT SOURCES)
    message(FATAL_ERROR
      "You must specify unit test source files using SOURCES")
  endif ()

  # Include the correct main file
  if (POLICY STREQUAL "MPI")
    set(_main_file ${CMAKE_SOURCE_DIRECTORY}/cmake/unit_main_mpi.cc)
  else ()
    set(_main_file ${CMAKE_SOURCE_DIRECTORY}/cmake/unit_main_serial.cc)
  endif ()

  list(INSERT ${SOURCES} 0 ${_main_file}) 
  
  # Set up the unit test executable and dependencies
  
  add_executable(${name} ${SOURCES})
  target_link_libraries(${name} PRIVATE ${GTEST_LIBRARIES})
  target_compile_definitions(${name} PRIVATE ${DEFINES})
  
  if (LIBRARIES)
    target_link_libraries(${name} PRIVATE ${LIBRARIES})
  endif ()
  
  if (INPUTS)
    set(_INPUT_FILES)
    foreach (input ${INPUTS})
      get_filename_component(_OUTPUT_NAME ${input} NAME)
      get_filename_component(_PATH ${input} ABSOLUTE)
      configure_file(${_PATH} ${PROJECT_BINARY_DIR}/${_OUTPUT_NAME})
      list(APPEND _INPUT_FILES ${PROJECT_BINARY_DIR}/${_OUTPUT_NAME})
    endforeach ()

    add_custom_target(${name}_inputs DEPENDS ${_INPUT_FILES})
    add_dependencies(${name} ${name}_inputs)
  endif ()

  
  if (POLICY STREQUAL "MPI")
    if (NOT THREADS)
      set(THREADS 1)
    endif ()
    
    foreach (instance ${THREADS})
      if (ENABLE_JENKINS_OUTPUT)
        set(_OUTPUT ${name}_${instance}.xml)
        set(_GTEST_FLAGS "--gtest_output=xml:${_OUTPUT}")
      endif()
      
      add_test(NAME ${name}
        COMMAND
	${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${THREADS}
        $<TARGET_FILE:${name}> ${ARGUMENTS}
        ${GTEST_FLAGS}
        WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
    endforeach ()

  else()

    if (ENABLE_JENKINS_OUTPUT)
      set(_OUTPUT ${PROJECT_BINARY_DIR}/${name}.xml)
      set(GTEST_FLAGS "--gtest_output=xml:${_OUTPUT}")
    endif()

    add_test(
      NAME ${name}
      COMMAND
      $<TARGET_FILE:${name}>
      ${GTEST_FLAGS}
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
  endif ()

endfunction (wonton_add_unit)
