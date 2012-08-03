#------------------------------------------------------------------------------#
# libdetran/src/drivers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(DRIVERS_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET( DRIVERS_SRC 
  Execute.cc
  StupidParser.cc
)

#-------------------------------------#
# EXECUTABLE 
#-------------------------------------#

ADD_EXECUTABLE( detran 
                detran.cc 
)
TARGET_LINK_LIBRARIES( detran
                       drivers 
                       solvers
                       transport
                       materials
                       angle
                       geometry
                       utilities
                       ioutils
                       ${GPERFTOOLS_LIBRARIES}                        
)
INSTALL_TARGETS(/bin detran)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY( drivers 
             SHARED
             ${DRIVERS_SRC} 
)
TARGET_LINK_LIBRARIES( drivers 
                       solvers
                       transport
                       materials
                       angle
                       geometry
                       utilities
)
INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib drivers)



#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#SET_SOURCE_FILES_PROPERTIES(detran_utilities.i PROPERTIES CPLUSPLUS ON)
#SWIG_ADD_MODULE(detran_utilities python detran_utilities.i )

#SWIG_LINK_LIBRARIES(detran_utilities ${PYTHON_LIBRARIES} utilities)

#INSTALL(FILES ${CMAKE_BINARY_DIR}/utilities/detran_utilities.py
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
#INSTALL(FILES ${CMAKE_BINARY_DIR}/utilities/_detran_utilities.so
#        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${DRIVERS_DIR}/test)
    add_subdirectory(${DRIVERS_DIR}/test)
endif()

