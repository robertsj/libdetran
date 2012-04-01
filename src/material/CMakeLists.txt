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
)

#-------------------------------------#
# PYTHON INTERFACE
#-------------------------------------#

IF( DETRAN_ENABLE_PYTHON )

SET_SOURCE_FILES_PROPERTIES(detran_materials.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(detran_materials python detran_materials.i )
SWIG_LINK_LIBRARIES(detran_materials
                    ${PYTHON_LIBRARIES} 
                    utilities 
                    materials)

ENDIF()

#-------------------------------------#
# TESTING
#-------------------------------------#

if( DETRAN_ENABLE_TEST )
    include_directories(${MATERIAL_DIR}/test)
    add_subdirectory(${MATERIAL_DIR}/test)
endif()

