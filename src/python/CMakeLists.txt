#------------------------------------------------------------------------------#
# src/python/CMakeLists.txt
#------------------------------------------------------------------------------#

# Install the helper initialization file
INSTALL_FILES(/python/detran __init__.py)

# add the utilities
add_subdirectory(pydetranutils)

