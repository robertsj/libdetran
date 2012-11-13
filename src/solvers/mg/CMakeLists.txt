# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MG_SRC
  ${SRC_DIR}/MGSolver.cc
  ${SRC_DIR}/MGTransportSolver.cc
  ${SRC_DIR}/MGSolverGS.cc
  ${SRC_DIR}/MGSolverGMRES.cc
  ${SRC_DIR}/MGDiffusionSolver.cc
  ${SRC_DIR}/MGTransportOperator.cc
  ${SRC_DIR}/DiffusionGainOperator.cc
  ${SRC_DIR}/DiffusionLossOperator.cc
  ${SRC_DIR}/MGPreconditioner.cc
  ${SRC_DIR}/MGDSA.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/solver/mg .hh)

