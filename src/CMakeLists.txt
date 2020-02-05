#------------------------------------------------------------------------------#
# libdetran/src/CMakeLists.txt
#------------------------------------------------------------------------------#

project(detran)
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)
set(DETRAN_DIR ${CMAKE_CURRENT_SOURCE_DIR})
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

#------------------------------------------------------------------------------#
# OPTIONS
#------------------------------------------------------------------------------#

# Set the default build type (debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug CACHE STRING
      "Choose the type of build, options are: Debug Release."
      FORCE)
endif()

# Options for shared library. On or Off.
set(BUILD_SHARED_LIBS YES CACHE BOOL        "Enable shared libraries.")
set(LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()

# Options for unit tests. On or Off.
set(DETRAN_ENABLE_TEST NO CACHE BOOL        "Enable the units tests.")

# Options for Python bindings. On or Off.
set(DETRAN_ENABLE_PYTHON NO CACHE BOOL      "Enable Python bindings via SWIG.")
# Options for OpenMP directives. On or Off.
set(DETRAN_ENABLE_OPENMP NO CACHE BOOL      "Enable OpenMP for multithreaded solves.")
# Options for PETSc. On or Off.
set(DETRAN_ENABLE_PETSC NO CACHE BOOL       "Build callow with PETSc.")
# Options for SLEPc. On or Off.
set(DETRAN_ENABLE_SLEPC NO CACHE BOOL       "Build callow SLEPc.")
# Options for gperftools. On or Off.
set(DETRAN_ENABLE_GPERFTOOLS NO CACHE BOOL  "Enable gperftools for profiling.")
# Options for silo. On or Off.
set(DETRAN_ENABLE_SILO NO CACHE BOOL        "Enable SILO.")
# Options for HDF5. On or Off.
set(DETRAN_ENABLE_HDF5 NO CACHE BOOL        "Enable HDF5.")
# Options for Boost. On or Off.
set(DETRAN_ENABLE_BOOST NO CACHE BOOL       "Enable BOOST.")
# Enable Fortran for Callow features
set(DETRAN_ENABLE_FORTRAN YES CACHE BOOL    "Enable Fortran.")

#------------------------------------------------------------------------------#
# CONFIGURATION
#------------------------------------------------------------------------------#

include(GetTheDate)
get_the_date(DETRAN_COMPILED_D DETRAN_COMPILED_M DETRAN_COMPILED_Y)

set(DETRAN_WINDOWS NO)
if(WIN32)
  set(DETRAN_WINDOWS YES)
endif()

set(DETRAN_WC_REVISION_HASH "unknown")
find_package(Git)
if(GIT_FOUND)
  GIT_WC_INFO(${CMAKE_CURRENT_SOURCE_DIR} "DETRAN")
endif()

# Write config file, install it, and include that directory for all source
set(DETRAN_CONFIG_IN  
    "${CMAKE_CURRENT_SOURCE_DIR}/detran_config/detran_config.hh.in")
set(DETRAN_CONFIG_OUT 
    "${PROJECT_BINARY_DIR}/detran_config/detran_config.hh")
configure_file(${DETRAN_CONFIG_IN} ${DETRAN_CONFIG_OUT})
include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_BINARY_DIR}/detran_config")
install(FILES       ${PROJECT_BINARY_DIR}/detran_config/detran_config.hh
        DESTINATION include)
install(FILES       "${CMAKE_CURRENT_SOURCE_DIR}/../README"
        DESTINATION .)
install(FILES 	        "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE"
      	DESTINATION .)

# Enable Fortran (on by default)
if(DETRAN_ENABLE_FORTRAN)
  enable_language (Fortran)
endif()

#------------------------------------------------------------------------------#
# EXTERNAL LIBRARIES
#------------------------------------------------------------------------------#

if(DETRAN_ENABLE_PYTHON)
  if(BUILD_SHARED_LIBS)
    find_package(SWIG REQUIRED)
    include(${SWIG_USE_FILE})
    find_package(PythonLibs REQUIRED)
    # find numpy/arrayobject.h
    exec_program(python 
                 ARGS            -c "'import numpy; print numpy.get_include()'"
                 OUTPUT_VARIABLE NPARRAYPATHHINT)
    find_path(NPARRAYPATH NAMES numpy/arrayobject.h PATHS ${NPARRAYPATHHINT})
    if (${NPARRAYPATH-NOTFOUND})
      message("!! Numpy could not be found and is needed for Python builds." 
              FATAL_ERROR)
    else()
      message("-- Found Numpy include: ${NPARRAYPATH}")
      include_directories(${NPARRAYPATH})
    endif()
    include_directories(${PYTHON_INCLUDE_PATH})
    include_directories(${PROJECT_BINARY_DIR}/detran_config)
    set(CMAKE_SWIG_FLAGS "-v")
  else()
    set(DETRAN_ENABLE_PYTHON OFF)
    message(STATUS "** Disabling Python since it required shared libraries.")
  endif()
endif()

if(DETRAN_ENABLE_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS ${OpenMP_CXX_FLAGS})
endif()

if(DETRAN_ENABLE_PETSC)
  find_package(PETSc REQUIRED)
  include_directories(${PETSC_INCLUDES})
  add_definitions(${PETSC_DEFINITIONS})
endif()

if(DETRAN_ENABLE_SLEPC)
  find_package(SLEPc REQUIRED)
  include_directories(${SLEPC_INCLUDES})
endif()

if(DETRAN_ENABLE_SILO)
  find_package(Silo)
  if(Silo_FOUND)
    include_directories(${Silo_INCLUDE_DIR})
  else()
    set(DETRAN_ENABLE_SILO OFF)
    message(STATUS "Disabling Silo since it was not found.")
  endif()
endif()

if(DETRAN_ENABLE_HDF5)
  find_package(HDF5)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    message("HDF5_INCLUDE_DIR = ${HDF5_INCLUDE_DIR}")
  else()
    set(DETRAN_ENABLE_HDF5 OFF)
    message(STATUS "Disabling HDF5 since it was not found.")
    # Need to reset them, since FindHDF5 fills them with "not founds".
    # which leads to errors downstream.
    set(HDF5_LIBRARIES "")
    set(HDF5_INCLUDE_DIR "")
  endif()
endif()

# This is a hack for now.  Probably worth generating a CMake finder thing.
if(DETRAN_ENABLE_GPERFTOOLS)
  set(GPERFTOOLS_LIBRARIES "/home/robertsj/opt/gperftools/lib/libprofiler.so")
  set(GPERFTOOLS_INCLUDE   "/home/robertsj/opt/gperftools/include")
  include_directories(${GPERFTOOLS_INCLUDE})
endif()

if(DETRAN_ENABLE_BOOST)
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREAD OFF)
  set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0")
  find_package(Boost 1.50 COMPONENTS serialization)
  if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
  else(Boost_FOUND)
    set(DETRAN_ENABLE_BOOST OFF)
    message(STATUS "Disabling BOOST since it was not found.")
    # Reset to avoid downstream errors due to NOTFOUND's.
    set(Boost_LIBRARIES "")
    set(Boost_INCLUDE_DIRS "")
  endif(Boost_FOUND)
endif()


#------------------------------------------------------------------------------#
# TESTING
#------------------------------------------------------------------------------#

OPTION (USE_TESTING "Build testing" ON)
if (DETRAN_ENABLE_TEST)
  set(MEMCHECK_COMMAND "/usr/bin/valgrind --leak-check=full")
  set(CTEST_INITIAL_CACHE "COVERAGE_COMMAND:FILEPATH=/usr/bin/gcov")
  enable_testing()
  include(CTest)
  # is there a better way?
  include_directories(utilities/test)
  include_directories(material/test)
  include_directories(geometry/test)
  include_directories(angle/test)
  include_directories(callow/test)
  include_directories(external_source/test)
  include_directories(discretization/test)
  include_directories(boundary/test)
  include_directories(transport/test)
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
include_directories(callow)
include_directories(external_source)
include_directories(orthog)
# level 2
include_directories(boundary)
include_directories(transport)
include_directories(kinetics)
# level 4
include_directories(solvers)
# level 5
include_directories(drivers)
#
include_directories(ioutils)
include_directories(postprocess)

# Set 
if(NOT CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR} CACHE STRING
      "Installation prefix; default is build directory."
      FORCE)
endif()

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include(SetupDetranLibrary)
add_subdirectory(utilities)
add_subdirectory(callow)
add_subdirectory(material)
add_subdirectory(angle)
add_subdirectory(geometry)
add_subdirectory(external_source)
add_subdirectory(orthog)
add_subdirectory(boundary)
add_subdirectory(transport)
add_subdirectory(ioutils)
add_subdirectory(kinetics)
add_subdirectory(solvers)
add_subdirectory(postprocess)
add_subdirectory(drivers)
add_subdirectory(cmake)

# Python things for installation
if(DETRAN_ENABLE_PYTHON)
  add_subdirectory(python)
endif()

#------------------------------------------------------------------------------#
# CPACK -- Packs up the build for installation elsewhere
#------------------------------------------------------------------------------#

if (PACKUP_DETRAN)
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Detran: a DETerministic TRANsport tool")
  set(CPACK_PACKAGE_NAME "Detran")
  set(CPACK_PACKAGE_VENDOR "Jeremy Roberts")
  set(CPACK_PACKAGE_VERSION "0.2.1")
  set(CPACK_PACKAGE_VERSION_MAJOR "0")
  set(CPACK_PACKAGE_VERSION_MINOR "2")
  set(CPACK_PACKAGE_VERSION_PATCH "1")
  set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE")
  set(CPACK_RESOURCE_FILE_README  "${CMAKE_CURRENT_SOURCE_DIR}/../README")
  if (WIN32)
    set(CPACK_NSIS_DISPLAY_NAME "Detran")
    set(CPACK_NSIS_HELP_LINK 
        "http:\\\\\\\\www.github.com/robertsj/libdetran")
    set(CPACK_NSIS_URL_INFO_ABOUT 
        "http:\\\\\\\\www.github.com/robertsj/libdetran")
    set(CPACK_NSIS_CONTACT "robertsj@mit.edu")
    set(CPACK_NSIS_MODIFY_PATH ON)
    set(CPACK_NSIS_EXTRA_INSTALL_COMMANDS 
        " \${EnvVarUpdate} \"$0\" \"LALALALALA\"  \"A\" \"HKCU\" \"${CPACK_PACKAGE_INSTALL_DIRECTORY}\\\\python\" ")
    set(CPACK_NSIS_MENU_LINKS "detran.bat;Detran")
  #elseif(APPLE)
  #  set(CPACK_GENERATOR "Bundle")
  #  include(AppleBundleMacro)
  #  ConfigureMacOSXBundlePlist("Detran", 
  #                             "_debug", 
  #                             "${CMAKE_CURRENT_SOURCE_DIR}/cmake/bundle.png"
  #                             "0.2.1")
  else()
    set(CPACK_GENERATOR "TGZ")
  endif()
  include(CPack)
endif()

#------------------------------------------------------------------------------#
# OUTPUT
#------------------------------------------------------------------------------#

message("---------------------------------------------------------------------")
message("              SUMMARY OF CMAKE CONFIGURATION")
message("---------------------------------------------------------------------")
message("++++ Date:        ${DETRAN_COMPILED_M}/${DETRAN_COMPILED_D}/${DETRAN_COMPILED_Y}")
message("++++ System:      ${CMAKE_SYSTEM} ${CMAKE_SYSTEM_PROCESSOR}")
message("++++ CXX:         ${CMAKE_CXX_COMPILER}")
message("++++ Build type:  ${CMAKE_BUILD_TYPE}")
message("++++ Install:     ${CMAKE_INSTALL_PREFIX}")

message("---------------------------------------------------------------------")
if(BUILD_SHARED_LIBS)
message("++++ Shared:      enabled")
else()
message("++++ Shared:      disabled")
endif()
if(DETRAN_ENABLE_PYTHON)
message("++++ Python:      enabled")
else()
message("++++ Python:      disabled")
endif()
if(DETRAN_ENABLE_OPENMP)
message("++++ OpenMP:      enabled")
else()
message("++++ OpenMP:      disabled")
endif()
if(DETRAN_ENABLE_PETSC)
message("++++ PETSc:       enabled")
else()
message("++++ PETSc:       disabled")
endif()
if(DETRAN_ENABLE_SLEPC)
message("++++ SLEPc:       enabled")
else()
message("++++ SLEPc:       disabled")
endif()
if(DETRAN_ENABLE_GPERFTOOLS)
message("++++ Profiling:   enabled")
else()
message("++++ Profiling:   disabled")
endif()
if(DETRAN_ENABLE_SILO)
message("++++ Silo:        enabled")
else()
message("++++ Silo:        disabled")
endif()
if(DETRAN_ENABLE_HDF5)
message("++++ HDF5:        enabled")
else()
message("++++ HDF5:        disabled")
endif()
if(DETRAN_ENABLE_BOOST)
message("++++ BOOST:       enabled")
else()
message("++++ BOOST:       disabled")
endif()
if(DETRAN_ENABLE_TEST)
message("++++ Testing:     enabled")
else()
message("++++ Testing:     disabled")
endif()
message("---------------------------------------------------------------------")
message("To use, remember to add:")
message("    ${CMAKE_INSTALL_PREFIX}/bin    to PATH")
message("    ${CMAKE_INSTALL_PREFIX}/lib    to LD_LIBRARY_PATH")
message("    ${CMAKE_INSTALL_PREFIX}/python to PYTHONPATH")
message("---------------------------------------------------------------------")


