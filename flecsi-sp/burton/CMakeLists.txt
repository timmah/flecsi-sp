#~----------------------------------------------------------------------------~#
# Copyright (c) 2016 Los Alamos National Laboratory, LLC
# All rights reserved
#~----------------------------------------------------------------------------~#

file(GLOB headers RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.h
  ${CMAKE_CURRENT_SOURCE_DIR}/detail/*.h
)

file(GLOB sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/*.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/detail/*.cc
)

set(burton_HEADERS
  ${headers}
  PARENT_SCOPE
)


set(burton_SOURCES
  ${sources} 
  PARENT_SCOPE
)


# get a relative path
file(RELATIVE_PATH DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}" ${FleCSI_SP_DATA_DIR})
