#------------------------------------------------------------------------------#
# src/diffusion/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    BaseOperator.cc
    OneGroupLossOperator.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(diffusion
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( diffusion
                       materials
                       geometry
                       utilities
                       ${PETSC_LIBRARIES}
                       ${GPERFTOOLS_LIBRARIES}                       
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib diffusion)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

