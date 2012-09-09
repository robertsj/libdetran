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

ADD_EXECUTABLE( xdetran 
                detran.cc 
)
TARGET_LINK_LIBRARIES( xdetran
                       detran 
                       solvers
                       transport
                       material
                       angle
                       geometry
                       utilities
                       ioutils
                       ${GPERFTOOLS_LIBRARIES}                        
)
INSTALL_TARGETS(/bin xdetran)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

# This is the detran library to be linked against.
ADD_LIBRARY( detran 
             SHARED
             ${DRIVERS_SRC} 
)
TARGET_LINK_LIBRARIES( detran 
                       solvers
                       transport
                       material
                       angle
                       geometry
                       utilities
                       ioutils
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib detran)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${DRIVERS_DIR}/test)
    add_subdirectory(${DRIVERS_DIR}/test)
endif()

