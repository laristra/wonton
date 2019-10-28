#[[
This file is part of the Ristra wonton project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]


# Add header files

if (FleCSI_DIR)
  add_library(flecsi_state_wrapper INTERFACE)
  
  set(flecsi_state_wrapper_HEADERS
    flecsi/flecsi_state_wrapper.h
    )

  target_sources(flecsi_state_wrapper INTERFACE ${flecsi_state_wrapper_HEADERS})

  # Add flecsi dependencies
  target_link_libraries(mesh INTERFACE ${FleCSI_LIBRARIES} ${FleCSISP_LIBRARIES})

  # Have to add FleCSI state wrapper tests
endif (FleCSI_DIR)
