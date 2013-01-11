# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(WG_SRC
  ${SRC_DIR}/WGSolver.cc
  ${SRC_DIR}/WGSolverSI.cc
  ${SRC_DIR}/WGSolverGMRES.cc
  ${SRC_DIR}/WGTransportOperator.cc
  ${SRC_DIR}/WGDiffusionLossOperator.cc
  ${SRC_DIR}/WGPreconditioner.cc
  ${SRC_DIR}/PC_DSA.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/solvers/wg .hh)

