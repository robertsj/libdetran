# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(MG_SRC
  ${SRC_DIR}/MGSolver.cc
  ${SRC_DIR}/MGTransportSolver.cc
  ${SRC_DIR}/MGSolverGS.cc
  ${SRC_DIR}/MGSolverGMRES.cc
  ${SRC_DIR}/MGSolverCMFD.cc
  ${SRC_DIR}/MGDiffusionSolver.cc
  ${SRC_DIR}/MGTransportOperator.cc
  ${SRC_DIR}/DiffusionGainOperator.cc
  ${SRC_DIR}/DiffusionLossOperator.cc
  ${SRC_DIR}/MGPreconditioner.cc
  ${SRC_DIR}/MGCoarseMeshPreconditioner.cc
  ${SRC_DIR}/MGDSA.cc
  ${SRC_DIR}/MGCMDSA.cc
  ${SRC_DIR}/MGTCDSA.cc
  ${SRC_DIR}/MGScatterFissionOperator.cc
  ${SRC_DIR}/MGSweepOperator.cc
  ${SRC_DIR}/SpectrumGS.cc
  ${SRC_DIR}/SpectrumFS.cc
  ${SRC_DIR}/SpectrumPinCell.cc
  ${SRC_DIR}/CMFDLossOperator.cc
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/solvers/mg .hh)

