# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(PRECONDITIONER_SRC
  ${SRC_DIR}/Preconditioner.cc
  ${SRC_DIR}/PCJacobi.cc
  ${SRC_DIR}/PCILU0.cc
  ${SRC_DIR}/PCShell.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/preconditioner .hh)


