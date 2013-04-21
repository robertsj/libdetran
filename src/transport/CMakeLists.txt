#------------------------------------------------------------------------------#
# src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

set(LIB_NAME "transport")

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

set(SRC
    BoundaryTally.cc
    CoarseMesh.cc
    CurrentTally.cc
    FissionSource.cc
    Homogenize.cc
    ScatterSource.cc
    State.cc
    Sweeper.cc
    Sweeper1D.cc
    Sweeper2D.cc
    Sweeper3D.cc
    Sweeper2DMOC.cc
    # discretization
    Equation_DD_1D.cc
    Equation_DD_2D.cc
    Equation_DD_3D.cc
    Equation_SC_1D.cc
    Equation_SC_2D.cc
    Equation_SD_1D.cc
    Equation_SD_2D.cc
    Equation_SC_MOC.cc
)

#-----------------------------------------------------------------------------#
# LINKED LIBRARIES
#-----------------------------------------------------------------------------#

set(LINKED_LIBS
    angle
    boundary
    external_source
    geometry
    material
    utilities
)

#-----------------------------------------------------------------------------#
# LIBRARY
#-----------------------------------------------------------------------------#

setup_detran_library(${LIB_NAME} "${SRC}" "${LINKED_LIBS}")
