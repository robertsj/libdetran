#------------------------------------------------------------------------------#
# libdetran/src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(TRANSPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(TRANSPORT_SRC
    State.cc
    Equation_DD_2D.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(transport
            SHARED
            ${TRANSPORT_SRC} 
)
TARGET_LINK_LIBRARIES( transport
                       angle
                       materials
                       geometry
                       utilities
)

ADD_EXECUTABLE(test_State test_State.cc)
TARGET_LINK_LIBRARIES( test_State
                       transport
                       geometry
                       angle
                       utilities )

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