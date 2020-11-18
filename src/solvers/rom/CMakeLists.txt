# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(ROM_SRC
  ${SRC_DIR}/ROMBasis.cc
  ${SRC_DIR}/OperatorProjection.cc
  ${SRC_DIR}/ROMSolver.cc
  ${SRC_DIR}/KineticMatrices.cc
  ${SRC_DIR}/TransientSolver.cc
  PARENT_SCOPE
)
