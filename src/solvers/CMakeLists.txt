#------------------------------------------------------------------------------#
# libdetran/src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SOLVERS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SOLVERS_SRC
	Eigensolver.cc
	EigenSLEPc.cc
    InnerIteration.cc
    SourceIteration.cc
    InnerGMRES.cc
    GaussSeidel.cc
    KrylovMG.cc
    PowerIteration.cc
    PreconditionerWG.cc
    PreconditionerMG.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(solvers
            SHARED
            ${SOLVERS_SRC} 
)
TARGET_LINK_LIBRARIES( solvers
                       diffusion
                       transport
                       angle
                       materials
                       geometry
                       utilities
                       ${PETSC_LIBRARIES}
                       ${SLEPC_LIBRARIES}  
                       ${GPERFTOOLS_LIBRARIES}
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib solvers)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SOLVERS_DIR}/test)
    add_subdirectory(${SOLVERS_DIR}/test)
endif()
