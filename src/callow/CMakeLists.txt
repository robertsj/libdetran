#------------------------------------------------------------------------------#
# callow - a naive sparse linear algebra library
#------------------------------------------------------------------------------#

SET(LIB_NAME "callow")

SET(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories(${SRC_DIR})
include_directories(utils)
include_directories(vector)
include_directories(matrix)
include_directories(solver)
include_directories(preconditioner)
add_subdirectory(utils)
add_subdirectory(vector)
add_subdirectory(matrix)
add_subdirectory(solver)
add_subdirectory(preconditioner)

#------------------------------------------------------------------------------#
# LIBRARY
#------------------------------------------------------------------------------#

add_library(callow
            SHARED
            ${UTILS_SRC}
            ${VECTOR_SRC}
            ${MATRIX_SRC}
            ${SOLVER_SRC}
            ${PRECONDITIONER_SRC}
)
target_link_libraries(callow 
                      utilities
                      ${PETSC_LIBRARIES}
                      ${SLEPC_LIBRARIES}  
                      ${GPERFTOOLS_LIBRARIES}
)
install_targets(/lib callow)
install_files(/include/callow .hh)


#-----------------------------------------------------------------------------#
# TESTING
#-----------------------------------------------------------------------------#

if(DETRAN_ENABLE_TEST)
    include_directories(${SRC_DIR}/test)
    add_subdirectory(${SRC_DIR}/test)
endif()

