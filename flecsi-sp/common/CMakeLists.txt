#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(common_HEADERS
  constants.h	
  exceptions.h
  parse_arguments.h
  specialization_init.h  
  types.h
  utils.h
  PARENT_SCOPE
)

set(common_SOURCES
  exceptions.cc
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
