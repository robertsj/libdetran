#------------------------------------------------------------------------------#
# src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

set(LIB_NAME "solvers")

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(wg)
include_directories(mg)
include_directories(eigen)
include_directories(time)
add_subdirectory(wg)
add_subdirectory(mg)
add_subdirectory(eigen)
add_subdirectory(time)
set(SRC
    FixedSourceManager.cc
    EigenvalueManager.cc
    SweepOperator.cc
    Solver.cc
    ${EIGEN_SRC}
    ${MG_SRC}
    ${WG_SRC}
    ${TIME_SRC}
)

#-----------------------------------------------------------------------------#
# LINKED LIBRARIES
#-----------------------------------------------------------------------------#

set(LINKED_LIBS
    callow
    transport
    angle
    material
    geometry
    utilities
    ioutils
    kinetics
)

#-----------------------------------------------------------------------------#
# LIBRARY
#-----------------------------------------------------------------------------#

setup_detran_library(${LIB_NAME} "${SRC}" "${LINKED_LIBS}")
