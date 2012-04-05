#------------------------------------------------------------------------------#
# libdetran/src/solvers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(TRANSPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(TRANSPORT_SRC
    InnerIteration.cc
    SourceIteration.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(solvers
            SHARED
            ${TRANSPORT_SRC} 
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

SET_SOURCE_FILES_PROPERTIES(detran_transport.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_transport python detran_transport.i )
SWIG_LINK_LIBRARIES(detran_transport
                    ${PYTHON_LIBRARIES} 
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
    include_directories(${TRANSPORT_DIR}/test)
    add_subdirectory(${TRANSPORT_DIR}/test)
endif()