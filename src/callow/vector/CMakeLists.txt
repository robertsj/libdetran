# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(VECTOR_SRC
  ${SRC_DIR}/Vector.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/vector .hh)


