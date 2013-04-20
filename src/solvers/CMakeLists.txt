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
include_directories(time)
add_subdirectory(wg)
add_subdirectory(mg)
add_subdirectory(eigen)
add_subdirectory(time)
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
            ${TIME_SRC}
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      callow
                      transport
                      angle
                      material
                      geometry
                      utilities
                      ioutils
                      kinetics
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
                    transport
                    
                    )

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/${LIB_NAME}.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/_${LIB_NAME}.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
        
endif()        
