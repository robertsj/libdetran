#------------------------------------------------------------------------------#
# src/kinetics/CMakeLists.txt
#------------------------------------------------------------------------------#

set(LIB_NAME "kinetics")

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

set(SRC
    KineticsMaterial.cc
    LinearMaterial.cc
    TimeDependentMaterial.cc
    PyTimeDependentMaterial.cc
    Precursors.cc
    SyntheticSource.cc
    SyntheticMomentSource.cc
    SyntheticDiscreteSource.cc
    LinearExternalSource.cc
    PulsedExternalSource.cc
    MultiPhysics.cc
)

#-----------------------------------------------------------------------------#
# LINKED LIBRARIES
#-----------------------------------------------------------------------------#

set(LINKED_LIBS
    material
    transport
)

#-----------------------------------------------------------------------------#
# LIBRARY
#-----------------------------------------------------------------------------#

setup_detran_library(${LIB_NAME} "${SRC}" "${LINKED_LIBS}")
