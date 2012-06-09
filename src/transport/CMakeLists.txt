#------------------------------------------------------------------------------#
# libdetran/src/transport/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(TRANSPORT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(TRANSPORT_SRC
    Acceleration.cc
    Boundary.cc
    CMR.cc
    Equation_DD_1D.cc
    Equation_DD_2D.cc
    Equation_DD_3D.cc
    Equation_SC_2D.cc
    Equation_SD_1D.cc
    Equation_SD_2D.cc
    FissionSource.cc
    ReactionRates.cc
    ScatterSource.cc
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
                       angle
                       materials
                       geometry
                       utilities
                       ${GPERFTOOLS_LIBRARIES}
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib transport)

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

#IF( DETRAN_ENABLE_PYTHON )

#SET_SOURCE_FILES_PROPERTIES(detran_transport.i PROPERTIES CPLUSPLUS ON)
#SWIG_ADD_MODULE(detran_transport python detran_transport.i )
#SWIG_LINK_LIBRARIES(detran_transport
#                    ${PYTHON_LIBRARIES} 
#                    transport 
#                    angle
#                    geometry
#                    materials
#                    utilities)

#INSTALL(FILES ${CMAKE_BINARY_DIR}/transport/detran_transport.py
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
#INSTALL(FILES ${CMAKE_BINARY_DIR}/transport/_detran_transport.so
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)

#ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${TRANSPORT_DIR}/test)
    add_subdirectory(${TRANSPORT_DIR}/test)
endif()
