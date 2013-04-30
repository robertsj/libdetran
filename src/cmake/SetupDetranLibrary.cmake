# - Function for automatic generation of a Detran library
#   given the name, sources, and linking libraries.
#
# The SETUP_DETRAN_LIBRARY function eliminates the boilerplate
# associated with adding a library.
#

include(GenerateExportHeaderModified)

function(SETUP_DETRAN_LIBRARY LIB_NAME LIB_SOURCE LIB_LINKS)

#message(" SETTING UP DETRAN LIBRARY = ${LIB_NAME}")
#message(" SOURCE = ${LIB_SOURCE}")
#message(" LINKS  = ${LIB_LINKS}")

# Add the target library and set the linked libraries
add_library(${LIB_NAME} ${LIB_TYPE} ${SRC})
target_link_libraries(${LIB_NAME} ${LIB_LINKS})

# Generate a header for Windows support
string(TOUPPER ${LIB_NAME} CAP_LIB_NAME)
generate_export_header(${LIB_NAME}
                       BASE_NAME         ${CAP_LIB_NAME}
                       EXPORT_MACRO_NAME ${CAP_LIB_NAME}_EXPORT
                       EXPORT_FILE_NAME  ${LIB_NAME}_export.hh
                       STATIC_DEFINE     ${CAP_LIB_NAME}_BUILT_AS_STATIC
)

# Install headers and other includes
FILE(GLOB_RECURSE files 
     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_SOURCE_DIR}/*.hh
     ${CMAKE_CURRENT_SOURCE_DIR}/*.i
)
install(FILES       ${files} 
        DESTINATION include/${LIB_NAME})
install(TARGETS     ${LIB_NAME} 
        DESTINATION lib)
install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/${LIB_NAME}_export.hh
        DESTINATION include/${LIB_NAME})


# Setup for testing
if(DETRAN_ENABLE_TEST)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/test)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/test)
endif()

# Setup for Python front end
if(DETRAN_ENABLE_PYTHON)
  # SWIG setup
  set_source_files_properties(detran_${LIB_NAME}.i PROPERTIES CPLUSPLUS ON)
  swig_add_module(${LIB_NAME} python detran_${LIB_NAME}.i)
  swig_link_libraries(${LIB_NAME} ${PYTHON_LIBRARIES} ${LIB_NAME})
  # Install the module and library
  install(FILES ${CMAKE_BINARY_DIR}/${LIB_NAME}/${LIB_NAME}.py
          DESTINATION python/detran)
  set(PYEXT ".so")
  set(PYPATH ${CMAKE_CURRENT_BINARY_DIR})
  if (MSVC)
    set(PYEXT ".pyd")
    set(PYPATH "${CMAKE_CURRENT_BINARY_DIR}/Release")
  endif()
  install(FILES       ${PYPATH}/_${LIB_NAME}${PYEXT}
          DESTINATION python/detran)
endif()

endfunction()
