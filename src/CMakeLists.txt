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

# Options for gperftools. On or Off.
SET( DETRAN_ENABLE_GPERFTOOLS NO CACHE BOOL
  "Enable Google Performance Tools for profiling."
)

# Options for silo. On or Off.
SET( DETRAN_ENABLE_SILO NO CACHE BOOL
  "Enable SILO."
)

# Options for HDF5. On or Off.
SET( DETRAN_ENABLE_HDF5 NO CACHE BOOL
  "Enable HDF5."
)

# Options for Boost. On or Off.
SET( DETRAN_ENABLE_BOOST NO CACHE BOOL
  "Enable BOOST."
)

# Write config file, install it, and include that directory for all source
set(DETRAN_CONFIG_IN  "${CMAKE_CURRENT_SOURCE_DIR}/config/detran_config.h.in")
set(DETRAN_CONFIG_OUT "${PROJECT_BINARY_DIR}/config/detran_config.h")
message("DETRAN_CONFIG_IN = ${DETRAN_CONFIG_IN}")
configure_file(${DETRAN_CONFIG_IN} ${DETRAN_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}/config")
install(FILES ${PROJECT_BINARY_DIR}/config/detran_config.h
        DESTINATION include)

#------------------------------------------------------------------------------#
# EXTERNAL LIBRARIES
#------------------------------------------------------------------------------#

IF( DETRAN_ENABLE_PYTHON )
  find_package(SWIG REQUIRED)
  INCLUDE(${SWIG_USE_FILE})
  find_package(PythonLibs REQUIRED)
  include_directories(${PYTHON_INCLUDE_PATH})
  include_directories(${PROJECT_BINARY_DIR}/config)
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

if( DETRAN_ENABLE_SILO)
  find_package(Silo)
  if(Silo_FOUND)
    include_directories(${Silo_INCLUDE_DIR})
  else()
      set( DETRAN_ENABLE_SILO OFF )
      MESSAGE( STATUS "Disabling Silo since it was not found." )
  endif()
endif()

if( DETRAN_ENABLE_HDF5)
  find_package(HDF5)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
  else()
      set( DETRAN_ENABLE_HDF5 OFF )
      MESSAGE( STATUS "Disabling HDF5 since it was not found." )
      # Need to reset them, since FindHDF5 fills them with "not founds", 
      # which leads to errors downstream.
      set(HDF5_LIBRARIES "")
      set(HDF_INCLUDE_DIR "")
  endif()
endif()

# This is a hack for now.  Probably worth generating a CMake finder thing.
IF( DETRAN_ENABLE_GPERFTOOLS )
  set(GPERFTOOLS_LIBRARIES "/home/robertsj/opt/gperftools/lib/libprofiler.so")
  set(GPERFTOOLS_INCLUDE   "/home/robertsj/opt/gperftools/include")
  include_directories(${GPERFTOOLS_INCLUDE})
ENDIF()

if( DETRAN_ENABLE_BOOST)
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREAD OFF)
  set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0")
  find_package(Boost 1.50 COMPONENTS serialization)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  else(Boost_FOUND)
    set( DETRAN_ENABLE_BOOST OFF )
    MESSAGE( STATUS "Disabling BOOST since it was not found." )
    # Need to reset them, since FindHDF5 fills them with "not founds", 
    # which leads to errors downstream.
    set(Boost_LIBRARIES "")
    set(Boost_INCLUDE_DIRS "")
  endif(Boost_FOUND)
endif()


#------------------------------------------------------------------------------#
# TESTING
#------------------------------------------------------------------------------#

OPTION ( USE_TESTING "Build testing" ON )
if ( DETRAN_ENABLE_TEST )
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
       /usr/bin/valgrind CACHE STRING 
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
    include_directories(external_source/test)
    include_directories(discretization/test)
    include_directories(boundary/test)
    include_directories(transport/test)
    include_directories(diffusion/test)
    include_directories(solvers/test)
    include_directories(drivers/test)
    include_directories(postprocess/test)
    include_directories(ioutils/test)
endif()


#------------------------------------------------------------------------------#
# INSTALL
#------------------------------------------------------------------------------#

include_directories(${DETRAN_DIR})
#
include_directories(python)
# level 0
include_directories(utilities)
# level 1
include_directories(angle)
include_directories(geometry)
include_directories(material)
include_directories(external_source)
# level 2
include_directories(discretization)
include_directories(boundary)
include_directories(transport)
# level 3
include_directories(diffusion)
# level 4
include_directories(solvers)
# level 5
include_directories(drivers)
#
include_directories(ioutils)
include_directories(postprocess)

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
add_subdirectory(external_source)
add_subdirectory(discretization)
add_subdirectory(boundary)
add_subdirectory(transport)
add_subdirectory(diffusion)
add_subdirectory(solvers)
add_subdirectory(drivers)
add_subdirectory(ioutils)
add_subdirectory(postprocess)

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("---------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION")
message("---------------------------------------------------------------")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")

message("---------------------------------------------------------------")
if( DETRAN_ENABLE_PYTHON )
message("++++ Python:      enabled")
else()
message("++++ Python:      disabled")
endif()
if( DETRAN_ENABLE_OPENMP )
message("++++ OpenMP:      enabled")
else()
message("++++ OpenMP:      disabled")
endif()
if( DETRAN_ENABLE_PETSC )
message("++++ PETSc:       enabled")
else()
message("++++ PETSc:       disabled")
endif()
if( DETRAN_ENABLE_SLEPC )
message("++++ SLEPc:       enabled")
else()
message("++++ SLEPc:       disabled")
endif()
if( DETRAN_ENABLE_GPERFTOOLS )
message("++++ Profiling:   enabled")
else()
message("++++ Profiling:   disabled")
endif()
if( DETRAN_ENABLE_SILO )
message("++++ Silo:        enabled")
else()
message("++++ Silo:        disabled")
endif()
if( DETRAN_ENABLE_HDF5 )
message("++++ HDF5:        enabled")
else()
message("++++ HDF5:        disabled")
endif()
if( DETRAN_ENABLE_BOOST )
message("++++ BOOST:       enabled")
else()
message("++++ BOOST:       disabled")
endif()
if( DETRAN_ENABLE_TEST )
message("++++ Testing:     enabled")
else()
message("++++ Testing:     disabled")
endif()
message("---------------------------------------------------------------")


