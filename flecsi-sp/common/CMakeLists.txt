#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(common_HEADERS
  types.h
  exception.h
  PARENT_SCOPE
)

set(common_SOURCES
  exception.cc
  PARENT_SCOPE
)

#cinch_install_headers(
#  FILES ${common_HEADERS}
#  DESTINATION include/flecsi-sp/common
#)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
