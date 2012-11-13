#------------------------------------------------------------------------------#
# src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

set(LIB_NAME "solvers")

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

# Include
include_directories(${SRC_DIR})
include_directories(wg)
include_directories(mg)
include_directories(eigen)
add_subdirectory(wg)
add_subdirectory(mg)
add_subdirectory(eigen)
set(SRC
    FixedSourceManager.cc
    EigenvalueManager.cc
    SweepOperator.cc
    Solver.cc
)

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
            ${WG_SRC}
            ${MG_SRC}
            ${EIGEN_SRC}
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      callow
                      diffusion
                      transport
                      angle
                      material
                      geometry
                      utilities
                      ioutils
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
