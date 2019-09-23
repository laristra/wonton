#[[
This file is part of the Ristra portage project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/wonton/blob/master/LICENSE
]]

add_executable(multidef_check main.cc dummy.cc)
target_link_libraries(multidef_check wonton ${LAPACKE_LIBRARIES})

