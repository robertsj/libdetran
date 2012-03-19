#------------------------------------------------------------------------------#
# libdetran/src/material/CMakeLists.txt
#------------------------------------------------------------------------------#

SET( MATERIAL_DIR ${CMAKE_CURRENT_SOURCE_DIR})

SET( MATERIAL_SRC
# driver
    driver/GlobalInput.cc
)


#message( "source files = ${ERME_SRC}" )

# Make the library
ADD_LIBRARY( material
             ${MATERIAL_SRC} 
)
TARGET_LINK_LIBRARIES(utililities
)

# Build the serment executable
if( SERMENT_BUILD_SLABTRAN )
    ADD_EXECUTABLE( sermentx driver/Serment.cc )
    TARGET_LINK_LIBRARIES( sermentx 
                           erme
                           diff2d
                           utils
                           ${PETSC_LIBRARIES} 
                           ${SLEPC_LIBRARIES}
                           ${SILO_LIBRARIES} 
    )
endif()

