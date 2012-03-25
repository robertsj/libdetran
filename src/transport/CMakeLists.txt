#------------------------------------------------------------------------------#
# libdetran/src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(TRANSPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(TRANSPORT_SRC
    State.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(transport
            SHARED
            ${TRANSPORT_SRC} 
)
TARGET_LINK_LIBRARIES( transport
                       materials
                       geometry
                       utilities
)

#ADD_EXECUTABLE(test_geometry test_geometry.cc)
#TARGET_LINK_LIBRARIES( test_geometry
#                       geometry
#                       utilities )

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

SET_SOURCE_FILES_PROPERTIES(detran_transport.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_transport python detran_transport.i )
SWIG_LINK_LIBRARIES(detran_transport
                    ${PYTHON_LIBRARIES} 
                    transport 
                    geometry
                    materials
                    utilities)

ENDIF()
