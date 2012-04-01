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

# BOOST (only header libraries, and really just shared_ptr)
# find_package(Boost 1.42.0  REQUIRED)
# include_directories(${Boost_INCLUDE_DIRS})

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

OPTION ( USE_TESTING "Build testing" ON )
if ( SERMENT_ENABLE_TEST )
    SET ( CTEST_INITIAL_CACHE "
      COVERAGE_COMMAND:FILEPATH=/usr/bin/gcov
    " )
    # need for coverage -- probably a cleaner way to set this up.
    SET( CMAKE_CXX_FLAGS_DEBUG        
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING 
         "Choose"
         FORCE
    )
    SET( CMAKE_CXX_FLAGS        
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    SET( CMAKE_C_FLAGS_DEBUG          
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    SET( CMAKE_EXE_LINKER_FLAGS 
         "-g -O0  -fprofile-arcs -ftest-coverage" CACHE STRING
         "Choose"
         FORCE
    )
    SET( CTEST_MEMORYCHECK_COMMAND 
       /home/robertsj/opt/valgrind/3.6.1/bin/valgrind CACHE STRING 
       "Set the memory checking command."
       FORCE
    )
    SET( CTEST_MEMORYCHECK_COMMAND_OPTIONS 
       "--tool=memcheck --leak-check=full" CACHE STRING 
       "Set the memory checking command options."
       FORCE
    ) 
    ENABLE_TESTING()
    INCLUDE(CTest)
    # is there a better way?
    include_directories(utilities/test)
    include_directories(material/test)
    include_directories(geometry/test)
    include_directories(angle/test)
    include_directories(transport/test)
endif()

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories(${DETRAN_DIR})
include_directories(utilities)
include_directories(material)
include_directories(geometry)
include_directories(angle)
include_directories(transport)
add_subdirectory(utilities)
add_subdirectory(material)
add_subdirectory(geometry)
add_subdirectory(angle)
add_subdirectory(transport)

#-------- OUTPUT ------------------------------------------------------#

message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
if( DETRAN_ENABLE_PYTHON )
message("++++ Python interface enabled")
endif()
if( DETRAN_ENABLE_TEST )
message("++++ CTest testing enabled")
endif()


