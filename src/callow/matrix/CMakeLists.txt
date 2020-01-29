# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MATRIX_SRC
  ${SRC_DIR}/MatrixBase.cc 
  ${SRC_DIR}/Matrix.cc
  ${SRC_DIR}/MatrixShell.cc
  ${SRC_DIR}/MatrixDense.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/matrix .hh)

