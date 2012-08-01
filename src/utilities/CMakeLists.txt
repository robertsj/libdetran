#------------------------------------------------------------------------------#
# libdetran/src/utilities/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(UTILITIES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET( UTILITIES_SRC 
  InputDB.cc
  IO_HDF5.cc
  GenException.cc
  Warning.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY( utilities 
             SHARED
             ${UTILITIES_SRC} 
)

TARGET_LINK_LIBRARIES( utilities
                       ${GPERFTOOLS_LIBRARIES}
                       ${HDF5_LIBRARIES} ) 

#-------------------------------------#
# INSTALL
#-------------------------------------#

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib utilities)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${UTILITIES_DIR}/test)
    add_subdirectory(${UTILITIES_DIR}/test)
endif()

