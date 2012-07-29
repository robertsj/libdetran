#------------------------------------------------------------------------------#
# libdetran/src/postprocess/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(SRC
  ReactionRates.cc
  SiloOutput.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(postprocess
            SHARED
            ${SRC} 
)
TARGET_LINK_LIBRARIES( postprocess
                       materials
                       geometry
                       utilities
                       transport
                       ${Silo_LIBRARIES}                  
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib postprocess)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

