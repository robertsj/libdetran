#------------------------------------------------------------------------------#
# src/callow/CMakeLists.txt
#------------------------------------------------------------------------------#

SET(LIB_NAME "callow")

#------------------------------------------------------------------------------#
# SOURCE
#------------------------------------------------------------------------------#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
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

set(SRC
    ${UTILS_SRC}
    ${VECTOR_SRC}
    ${MATRIX_SRC}
    ${SOLVER_SRC}
    ${PRECONDITIONER_SRC}
)

#-----------------------------------------------------------------------------#
# LINKED LIBRARIES
#-----------------------------------------------------------------------------#

set(LINKED_LIBS
    utilities
    ${PETSC_LIBRARIES}
    ${SLEPC_LIBRARIES}  
    ${GPERFTOOLS_LIBRARIES}
)

#-----------------------------------------------------------------------------#
# LIBRARY
#-----------------------------------------------------------------------------#

setup_detran_library(${LIB_NAME} "${SRC}" "${LINKED_LIBS}")
 