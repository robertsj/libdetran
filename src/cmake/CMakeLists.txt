# Install CMake finders to facilitate using Detran 

file(GLOB_RECURSE files 
     RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cmake
)

install(FILES       ${files} 
        DESTINATION tools/CMake)
