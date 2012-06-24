#------------------------------------------------------------------------------#
# libdetran/src/python/CMakeLists.txt
#------------------------------------------------------------------------------#

SET_SOURCE_FILES_PROPERTIES(detran.i PROPERTIES CPLUSPLUS ON)
SWIG_ADD_MODULE(pydetran python detran.i )
SWIG_LINK_LIBRARIES(pydetran
                    ${PYTHON_LIBRARIES} 
                    #
                    utilities 
                    #
                    angle
                    geometry
                    materials
                    moc
                    #
                    transport
                    #
                    solvers
                    #
                    ${PETSC_LIBRARIES})

# Install the module and library
INSTALL(FILES ${CMAKE_BINARY_DIR}/python/pydetran.py
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
INSTALL(FILES ${CMAKE_BINARY_DIR}/python/_pydetran.so
        DESTINATION ${CMAKE_INSTALL_PREFIX}/python/detran)
        
# Install the helper initialization file
INSTALL_FILES(/python/detran __init__.py)
