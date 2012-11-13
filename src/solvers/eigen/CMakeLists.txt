# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(EIGEN_SRC
  ${SRC_DIR}/Eigensolver.cc
  ${SRC_DIR}/EigenArnoldi.cc
  ${SRC_DIR}/EigenPI.cc
  ${SRC_DIR}/EigenDiffusion.cc
  ${SRC_DIR}/EnergyIndependentEigenOperator.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/solver/wg .hh)

