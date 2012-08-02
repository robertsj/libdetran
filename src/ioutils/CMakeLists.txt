#------------------------------------------------------------------------------#
# libdetran/src/ioutils/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
  IO_HDF5.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(ioutils
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( ioutils
                       materials
                       geometry
                       utilities
                       ${HDF5_LIBRARIES}                
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib ioutils)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

