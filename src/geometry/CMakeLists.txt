#------------------------------------------------------------------------------#
# libdetran/src/geometry/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(GEOMETRY_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(GEOMETRY_SRC
    Mesh.cc
    Mesh1D.cc
    Mesh2D.cc
    Mesh3D.cc
    PinCell.cc
    Assembly.cc
    Core.cc
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

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib geometry)

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

INSTALL(FILES ${CMAKE_BINARY_DIR}/geometry/detran_geometry.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/geometry/_detran_geometry.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${GEOMETRY_DIR}/test)
    add_subdirectory(${GEOMETRY_DIR}/test)
endif()
