#------------------------------------------------------------------------------#
# libdetran/src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SOLVERS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SOLVERS_SRC
    InnerIteration.cc
    SourceIteration.cc
    GaussSeidel.cc
    PowerIteration.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(solvers
            SHARED
            ${SOLVERS_SRC} 
)
TARGET_LINK_LIBRARIES( solvers
                       transport
                       angle
                       materials
                       geometry
                       utilities
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib solvers)

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

SET_SOURCE_FILES_PROPERTIES(detran_solvers.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_solvers python detran_solvers.i )
SWIG_LINK_LIBRARIES(detran_solvers
                    ${PYTHON_LIBRARIES} 
                    solvers
                    transport 
                    angle
                    geometry
                    materials
                    utilities)
                    
INSTALL(FILES ${CMAKE_BINARY_DIR}/solvers/detran_solvers.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/solvers/_detran_solvers.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)          

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SOLVERS_DIR}/test)
    add_subdirectory(${SOLVERS_DIR}/test)
endif()