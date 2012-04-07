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

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib utilities)

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET_SOURCE_FILES_PROPERTIES(detran_utilities.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_utilities python detran_utilities.i )

SWIG_LINK_LIBRARIES(detran_utilities ${PYTHON_LIBRARIES} utilities)

INSTALL(FILES ${CMAKE_BINARY_DIR}/utilities/detran_utilities.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/utilities/_detran_utilities.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${UTILITIES_DIR}/test)
    add_subdirectory(${UTILITIES_DIR}/test)
endif()

