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

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SOLVERS_DIR}/test)
    add_subdirectory(${SOLVERS_DIR}/test)
endif()