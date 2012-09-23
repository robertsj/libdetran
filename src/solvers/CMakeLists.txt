#------------------------------------------------------------------------------#
# src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "solvers")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

SET(SRC
    DiffusionFixedSourceSolver.cc
    DiffusionGainOperator.cc
    DiffusionLossOperator.cc
    Eigensolver.cc
	EigenSLEPc.cc
    InnerIteration.cc
    SourceIteration.cc
    InnerGMRES.cc
    GaussSeidelMG.cc
    KrylovMG.cc
    PowerIteration.cc
    PreconditionerWG.cc
    PreconditionerMG.cc
)

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      callow
                      diffusion
                      transport
                      angle
                      material
                      geometry
                      utilities
                      ${Callow_LIBRARIES}
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
