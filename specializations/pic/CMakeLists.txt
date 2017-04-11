#~----------------------------------------------------------------------------~#
# Copyright (c) 2014 Los Alamos National Security, LLC
# All rights reserved.
#~----------------------------------------------------------------------------~#

set(pic_HEADERS
  config.h
  entity_types.h
  mesh.h
  types.h
  particles.h
  boundary.h
  PARENT_SCOPE
)

set(pic_SOURCES
  PARENT_SCOPE
)

cinch_add_unit(pic
  SOURCES test/pic.cc
  LIBRARIES ${FleCSI_LIBRARY}
)

#----------------------------------------------------------------------------~-#
# Formatting options for vim.
# vim: set tabstop=2 shiftwidth=2 expandtab :
#----------------------------------------------------------------------------~-#
