#------------------------------------------------------------------------------#
# libdetran/src/utilities/CMakeLists.txt
#------------------------------------------------------------------------------#

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

ADD_EXECUTABLE(test_InputDB test_InputDB.cc)
TARGET_LINK_LIBRARIES( test_InputDB
                       utilities )


#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET_SOURCE_FILES_PROPERTIES(detran_utilities.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_utilities python detran_utilities.i )

SWIG_LINK_LIBRARIES(detran_utilities ${PYTHON_LIBRARIES} utilities)

ENDIF()
