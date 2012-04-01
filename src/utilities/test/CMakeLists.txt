#------------------------------------------------------------------------------#
# libdetran/src/utilities/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran utilities tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Testing test_Testing.cc)
TARGET_LINK_LIBRARIES(test_Testing utilities)


#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_Testing_pass test_Testing 1)
ADD_TEST( test_Testing_fail test_Testing 2)
