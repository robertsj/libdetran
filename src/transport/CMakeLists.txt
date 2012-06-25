#------------------------------------------------------------------------------#
# libdetran/src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(TRANSPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(TRANSPORT_SRC
    Acceleration.cc
    Boundary.cc
    CMR.cc
    Equation_DD_1D.cc
    Equation_DD_2D.cc
    Equation_DD_3D.cc
    Equation_SC_2D.cc
    Equation_SD_1D.cc
    Equation_SD_2D.cc
    Equation_SC_MOC.cc
    FissionSource.cc
    ReactionRates.cc
    ScatterSource.cc
    State.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(transport
            SHARED
            ${TRANSPORT_SRC} 
)
TARGET_LINK_LIBRARIES( transport
                       angle
                       materials
                       geometry
                       moc
                       utilities
                       ${GPERFTOOLS_LIBRARIES}
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib transport)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${TRANSPORT_DIR}/test)
    add_subdirectory(${TRANSPORT_DIR}/test)
endif()
