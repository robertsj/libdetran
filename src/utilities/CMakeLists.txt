#------------------------------------------------------------------------------#
# libdetran/src/utilities/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(UTILITIES_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET( UTILITIES_SRC 
  InputDB.cc
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
if( DETRAN_ENABLE_GPERFTOOLS )
TARGET_LINK_LIBRARIES( ${GPERFTOOLS_LIBRARIES}
)
endif()

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

