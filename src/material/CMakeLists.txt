#------------------------------------------------------------------------------#
# libdetran/src/material/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(MATERIAL_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#-------------------------------------#
# SOURCE
#-------------------------------------#

SET(MATERIAL_SRC
    Material.cc
)

#-------------------------------------#
# SHARED LIBRARY
#-------------------------------------#

ADD_LIBRARY(materials
            SHARED
            ${MATERIAL_SRC} 
)
TARGET_LINK_LIBRARIES( materials
                       utilities
                       ${Boost_LIBRARIES}                     
)

INSTALL_FILES(/include .hh)
INSTALL_TARGETS(/lib materials)

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${MATERIAL_DIR}/test)
    add_subdirectory(${MATERIAL_DIR}/test)
endif()

