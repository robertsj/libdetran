#------------------------------------------------------------------------------#
# libdetran/src/geometry/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(GEOMETRY_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(GEOMETRY_SRC
    Mesh.cc
    Mesh2D.cc
    Mesh3D.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(geometry
            SHARED
            ${GEOMETRY_SRC} 
)
TARGET_LINK_LIBRARIES( geometry
                       utilities
)

ADD_EXECUTABLE(test_geometry test_geometry.cc)
TARGET_LINK_LIBRARIES( test_geometry
                       geometry
                       utilities )

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

SET_SOURCE_FILES_PROPERTIES(detran_geometry.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_geometry python detran_geometry.i )
SWIG_LINK_LIBRARIES(detran_geometry
                    ${PYTHON_LIBRARIES} 
                    utilities 
                    geometry)

ENDIF()
