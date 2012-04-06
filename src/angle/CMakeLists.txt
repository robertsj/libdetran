#------------------------------------------------------------------------------#
# libdetran/src/angle/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(ANGLE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(ANGLE_SRC
    GaussLegendre.cc
    MomentToDiscrete.cc
    Quadrature.cc
    QuadrupleRange.cc
    SphericalHarmonics.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(angle
            SHARED
            ${ANGLE_SRC} 
)
TARGET_LINK_LIBRARIES( angle
                       utilities
)

#ADD_EXECUTABLE(test_geometry test_geometry.cc)
#TARGET_LINK_LIBRARIES( test_angle
#                       angle
#                       utilities )

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

SET_SOURCE_FILES_PROPERTIES(detran_angle.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_angle python detran_angle.i )
SWIG_LINK_LIBRARIES(detran_angle
                    ${PYTHON_LIBRARIES} 
                    utilities 
                    angle)

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${ANGLE_DIR}/test)
    add_subdirectory(${ANGLE_DIR}/test)
endif()