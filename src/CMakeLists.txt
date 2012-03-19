#------------------------------------------------------------------------------#
# libdetran/src/CMakeLists.txt
#------------------------------------------------------------------------------#

project(detran)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(DETRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if( NOT CMAKE_BUILD_TYPE)
  set( CMAKE_BUILD_TYPE Debug CACHE STRING
       "Choose the type of build, options are: Debug Release."
       FORCE)
endif()

# Options for Python bindings. On or Off.
SET( DETRAN_ENABLE_PYTHON NO CACHE BOOL
  "Enable Python bindings via SWIG."
)

# Write config file, install it, and include that directory for all source
set(DETRAN_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/detran_config.h.in")
set(DETRAN_CONFIG_OUT "${PROJECT_BINARY_DIR}/config/detran_config.h")
message("DETRAN_CONFIG_IN = ${DETRAN_CONFIG_IN}")
configure_file(${DETRAN_CONFIG_IN} ${DETRAN_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/config")

#------------------------------------------------------------------------------#
# EXTERNAL LIBRARIES
#------------------------------------------------------------------------------#

IF( DETRAN_ENABLE_PYTHON )
  find_package(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH})
  set(CMAKE_SWIG_FLAGS "")
ENDIF()


#------------------------------------------------------------------------------#
# TESTING
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

#add_subdirectory(material)
add_subdirectory(utilities)
#add_subdirectory(geometry)
#include_directories(material)
include_directories(utilities)

#-------- OUTPUT ------------------------------------------------------#

message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
if( DETRAN_ENABLE_PYTHON )
message("++++ Python interface enabled")
endif()



