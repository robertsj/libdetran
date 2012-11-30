#------------------------------------------------------------------------------#
# src/kinetics/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "kinetics")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

SET(SRC
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
    LRA.cc
)

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      material
                      transport
)

INSTALL_FILES(/include/${LIB_NAME} .hh)
INSTALL_TARGETS(/lib ${LIB_NAME})

#-----------------------------------------------------------------------------#
# TESTING
#-----------------------------------------------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()
