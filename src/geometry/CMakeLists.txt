#------------------------------------------------------------------------------#
# src/geometry/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
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
            ${SRC} 
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
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()
