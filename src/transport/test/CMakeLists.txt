#------------------------------------------------------------------------------#
# libdetran/src/transport/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran geometry tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_Equation_DD_2D test_Equation_DD_2D.cc)
TARGET_LINK_LIBRARIES(test_Equation_DD_2D transport geometry materials utilities)

ADD_EXECUTABLE(test_Sweeper2D test_Sweeper2D.cc)
TARGET_LINK_LIBRARIES(test_Sweeper2D transport geometry materials utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_Equation_DD_2D_basic test_Equation_DD_2D 0)
ADD_TEST(test_Sweeper2D_basic      test_Sweeper2D      0)