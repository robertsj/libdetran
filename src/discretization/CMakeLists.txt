#------------------------------------------------------------------------------#
# src/discretization/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "discretization")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

SET(SRC
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
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      angle
                      geometry
                      material
                      utilities
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
