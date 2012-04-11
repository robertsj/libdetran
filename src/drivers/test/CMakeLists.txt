#------------------------------------------------------------------------------#
# libdetran/src/drivers/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran drivers tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_StupidParser test_StupidParser.cc)
TARGET_LINK_LIBRARIES(test_StupidParser drivers utilities geometry materials transport)

INSTALL(FILES ${CMAKE_SOURCE_DIR}/drivers/test/test.inp
        DESTINATION ${CMAKE_BINARY_DIR}/drivers/test/)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

# StupidParser
#ADD_TEST( test_Testing_pass  test_Testing 0)