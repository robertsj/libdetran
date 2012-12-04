#------------------------------------------------------------------------------#
# src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "transport")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

SET(SRC
    BoundaryTally.cc
    CoarseMesh.cc
    CurrentTally.cc
    FissionSource.cc
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
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      angle
                      boundary
                      external_source
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

#-----------------------------------------------------------------------------#
# PYTHON
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_PYTHON)

SET_SOURCE_FILES_PROPERTIES(detran_${LIB_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(${LIB_NAME} python detran_${LIB_NAME}.i )
SWIG_LINK_LIBRARIES(${LIB_NAME}
                    ${PYTHON_LIBRARIES}
                    ${LIB_NAME} 
                    utilities
                    )

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/${LIB_NAME}.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/_${LIB_NAME}.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
        
endif()        
