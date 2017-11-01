#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(io_HEADERS
  io_exodus.h
  write_binary.h
  vtk.h
)

set(io_SOURCES
  vtk.cc
)

cinch_install_headers(
  FILES ${io_HEADERS}
  DESTINATION include/flecsi-sp/io
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
