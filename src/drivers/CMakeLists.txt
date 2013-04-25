#------------------------------------------------------------------------------#
# src/drivers/CMakeLists.txt
#------------------------------------------------------------------------------#

set(EXE_NAME "xdetran")

#-----------------------------------------------------------------------------#
# SOURCE
#-----------------------------------------------------------------------------#

set(SRC 
    detran.cc
    Execute.cc
    StupidParser.cc
)

#-----------------------------------------------------------------------------#
# EXECUTABLE
#-----------------------------------------------------------------------------#

add_executable(${EXE_NAME} 
               ${SRC}
)
target_link_libraries(${EXE_NAME}
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
# TESTING
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_TEST)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif()

