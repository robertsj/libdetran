#------------------------------------------------------------------------------#
# src/ioutils/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "ioutils")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
    IO_HDF5.cc
    SiloOutput.cc
    PPMPlotter.cc
)

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      material
                      geometry
                      utilities
                      transport
                      ${HDF5_LIBRARIES}       
                      ${Silo_LIBRARIES}          
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

#-----------------------------------------------------------------------------#
# PYTHON
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_PYTHON)

SET_SOURCE_FILES_PROPERTIES(detran_${LIB_NAME}.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(${LIB_NAME} python detran_${LIB_NAME}.i )
SWIG_LINK_LIBRARIES(${LIB_NAME}
                    ${PYTHON_LIBRARIES}
                    ${LIB_NAME} 
                    )

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/${LIB_NAME}.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/_${LIB_NAME}.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
        
endif()        
