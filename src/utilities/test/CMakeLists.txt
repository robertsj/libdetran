#------------------------------------------------------------------------------#
# libdetran/src/utilities/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran utilities tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Testing test_Testing.cc)
TARGET_LINK_LIBRARIES(test_Testing utilities)

ADD_EXECUTABLE(test_InputDB test_InputDB.cc)
TARGET_LINK_LIBRARIES(test_InputDB utilities)

ADD_EXECUTABLE(test_MathUtilities test_MathUtilities.cc)
TARGET_LINK_LIBRARIES(test_MathUtilities utilities)

ADD_EXECUTABLE(test_Point   test_Point.cc)
TARGET_LINK_LIBRARIES(test_Point utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

# Testing
ADD_TEST( test_Testing_pass  test_Testing 0)
ADD_TEST( test_Testing_fail  test_Testing 1)
# InputDB
ADD_TEST( test_InputDB_basic test_InputDB 0)
# MathUtilities
ADD_TEST( test_MathUtilities test_norm_L2 0)
ADD_TEST( test_MathUtilities test_norm_L1 1)
ADD_TEST( test_MathUtilities test_norm_Linf 2)
ADD_TEST( test_MathUtilities test_vec_scale 3)
ADD_TEST( test_MathUtilities test_norm_residual_L2 4)
ADD_TEST( test_MathUtilities test_norm_residual_L1 5)
ADD_TEST( test_MathUtilities test_norm_residual_Linf 6)
# Point
ADD_TEST( test_Point         test_Point   0)
