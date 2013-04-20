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
INSTALL_FILES(/include/${LIB_NAME} .i)
INSTALL_TARGETS(/lib ${LIB_NAME})

#-----------------------------------------------------------------------------#
# TESTING
#-----------------------------------------------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

#-----------------------------------------------------------------------------#
# PYTHON
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_PYTHON)

SET_SOURCE_FILES_PROPERTIES(detran_${LIB_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(${LIB_NAME} python detran_${LIB_NAME}.i )
SWIG_LINK_LIBRARIES(${LIB_NAME}
                    ${PYTHON_LIBRARIES}
                    ${LIB_NAME} 
                    utilities
                    angle
                    )

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/${LIB_NAME}.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/_${LIB_NAME}.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)

endif()        