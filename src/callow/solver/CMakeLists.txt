# Set source
SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
SET(SOLVER_SRC
  ${SRC_DIR}/LinearSolver.cc
  ${SRC_DIR}/Richardson.cc
  ${SRC_DIR}/Jacobi.cc
  ${SRC_DIR}/GaussSeidel.cc
  ${SRC_DIR}/GMRES.cc
  ${SRC_DIR}/PetscSolver.cc
  ${SRC_DIR}/LinearSolverCreator.cc
  ${SRC_DIR}/SlepcSolver.cc
  ${SRC_DIR}/EigenSolver.cc
  ${SRC_DIR}/PowerIteration.cc
  ${SRC_DIR}/EigenSolverCreator.cc
  ${SRC_DIR}/Davidson.cc
  ${SRC_DIR}/Eispack.cc
  ${SRC_DIR}/Eispack.f90
  PARENT_SCOPE
)

# Install
INSTALL_FILES(/include/callow/solver .hh)
