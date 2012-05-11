#------------------------------------------------------------------------------#
# libdetran/src/angle/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran angle tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

# Quadratures
ADD_EXECUTABLE(test_GaussLegendre  test_GaussLegendre.cc)
ADD_EXECUTABLE(test_QuadrupleRange test_QuadrupleRange.cc)
# Other
ADD_EXECUTABLE(test_MomentToDiscrete test_MomentToDiscrete.cc)

TARGET_LINK_LIBRARIES(test_GaussLegendre  angle utilities)
TARGET_LINK_LIBRARIES(test_QuadrupleRange angle utilities)
TARGET_LINK_LIBRARIES(test_MomentToDiscrete angle utilities)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_GaussLegendre_basic  test_GaussLegendre 0)
ADD_TEST(test_QuadrupleRange_basic test_QuadrupleRange 0)
ADD_TEST(test_MomentToDiscrete_1D  test_MomentToDiscrete 0)
ADD_TEST(test_MomentToDiscrete_2D  test_MomentToDiscrete 0)

