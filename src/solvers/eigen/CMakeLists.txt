# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(EIGEN_SRC
  ${SRC_DIR}/Eigensolver.cc
  ${SRC_DIR}/EigenArnoldi.cc
  ${SRC_DIR}/EigenPI.cc
  ${SRC_DIR}/EigenDiffusion.cc
  ${SRC_DIR}/EigenGD.cc
  ${SRC_DIR}/EigenCMFD.cc
  ${SRC_DIR}/EnergyIndependentEigenOperator.cc
  ${SRC_DIR}/EnergyDependentEigenLHS.cc
  PARENT_SCOPE
)