#------------------------------------------------------------------------------#
# src/drivers/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "detran")
SET(EXE_NAME "xdetran")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

SET(SRC 
    Execute.cc
    StupidParser.cc
)

#-----------------------------------------------------------------------------#
# EXECUTABLE
#-----------------------------------------------------------------------------#

ADD_EXECUTABLE(${EXE_NAME} 
               detran.cc 
)
TARGET_LINK_LIBRARIES(${EXE_NAME}
                      ${LIB_NAME} 
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

#-----------------------------------------------------------------------------#
# SHARED LIBRARY
#-----------------------------------------------------------------------------#

# This is the detran library to be linked against.
ADD_LIBRARY(${LIB_NAME}
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES(${LIB_NAME}
                      solvers
                      transport
                      material
                      angle
                      geometry
                      utilities
                      ioutils
                      postprocess
)

INSTALL_FILES(/include/drivers/${LIB_NAME} .hh)
INSTALL_TARGETS(/lib ${LIB_NAME})

#-----------------------------------------------------------------------------#
# TESTING
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_TEST)
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

