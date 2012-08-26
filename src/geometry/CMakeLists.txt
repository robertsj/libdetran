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
    TrackDB.cc
    Tracker.cc
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
                       angle
                       ${Boost_LIBRARIES}        
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib geometry)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${GEOMETRY_DIR}/test)
    add_subdirectory(${GEOMETRY_DIR}/test)
endif()
