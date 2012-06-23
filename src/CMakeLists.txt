#------------------------------------------------------------------------------#
# libdetran/src/CMakeLists.txt
#------------------------------------------------------------------------------#

project(detran)
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
set(DETRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
LIST( APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" )

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

# Options for OpenMP directives. On or Off.
SET( DETRAN_ENABLE_OPENMC NO CACHE BOOL
  "Enable OpenMP for multithreaded solves."
)

# Options for PETSc. On or Off.
SET( DETRAN_ENABLE_PETSC NO CACHE BOOL
  "Enable PETSc."
)

# Options for SLEPc. On or Off.
SET( DETRAN_ENABLE_SLEPC NO CACHE BOOL
  "Enable SLEPc."
)

# Options for OpenMP directives. On or Off.
SET( DETRAN_ENABLE_GPERFTOOLS NO CACHE BOOL
  "Enable Google Performance Tools for profiling."
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

IF( DETRAN_ENABLE_OPENMP )
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
ENDIF()

if( DETRAN_ENABLE_PETSC)
  find_package(PETSc REQUIRED)
  include_directories(${PETSC_INCLUDES})
  add_definitions(${PETSC_DEFINITIONS})
endif()

if( DETRAN_ENABLE_SLEPC)
  find_package(SLEPc REQUIRED)
  include_directories(${SLEPC_INCLUDES})
endif()

# This is a hack for now.  Probably worth generating a CMake finder thing.
IF( DETRAN_ENABLE_GPERFTOOLS )
  set(GPERFTOOLS_LIBRARIES "/home/robertsj/opt/gperftools/lib/libprofiler.so")
  set(GPERFTOOLS_INCLUDE   "/home/robertsj/opt/gperftools/include")
  include_directories(${GPERFTOOLS_INCLUDE})
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
    include_directories(moc/test)
    include_directories(transport/test)
    include_directories(solvers/test)
    include_directories(drivers/test)
endif()


#------------------------------------------------------------------------------#
# INSTALL
#------------------------------------------------------------------------------#

include_directories(${DETRAN_DIR})
#
include_directories(python)
#
include_directories(utilities)
#
include_directories(angle)
include_directories(geometry)
include_directories(material)
include_directories(moc)
#
include_directories(transport)
#
include_directories(solvers)
#
include_directories(drivers)

# Set 
if( NOT CMAKE_INSTALL_PREFIX)
  set( CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
       "Installation prefix; default is build directory."
       FORCE)
endif()

# Python things for installation
if( DETRAN_ENABLE_PYTHON )
  add_subdirectory(python)
endif()

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

add_subdirectory(utilities)
add_subdirectory(material)
add_subdirectory(geometry)
add_subdirectory(angle)
add_subdirectory(moc)
add_subdirectory(transport)
add_subdirectory(solvers)
add_subdirectory(drivers)


#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")
if( DETRAN_ENABLE_PYTHON )
message("++++ Python:      enabled")
elseif()
message("++++ Python:      disabled")
endif()
if( DETRAN_ENABLE_OPENMP )
message("++++ OpenMP:      enabled")
elseif()
message("++++ OpenMP:      disabled")
endif()
if( DETRAN_ENABLE_PETSC )
message("++++ PETSc:       enabled")
elseif()
message("++++ PETSc:       disabled")
endif()
if( DETRAN_ENABLE_SLEPC )
message("++++ SLEPc:       enabled")
elseif()
message("++++ SLEPc:       disabled")
endif()
if( DETRAN_ENABLE_GPERFTOOLS )
message("++++ Profiling:   enabled")
elseif()
message("++++ Profiling:   disabled")
endif()
if( DETRAN_ENABLE_TEST )
message("++++ Testing:     enabled")
elseif()
message("++++ Testing:     disabled")
endif()


