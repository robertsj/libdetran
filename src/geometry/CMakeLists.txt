#------------------------------------------------------------------------------#
# src/geometry/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "geometry")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

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

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      angle
                      utilities
)

INSTALL_FILES(/include/${LIB_NAME} .hh)
INSTALL_TARGETS(/lib ${LIB_NAME})

#-----------------------------------------------------------------------------#
# TESTING
#-----------------------------------------------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()
