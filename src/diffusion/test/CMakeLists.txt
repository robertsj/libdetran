#------------------------------------------------------------------------------#
# /src/diffusion/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran diffusion tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_GainOperator        test_GainOperator.cc)
TARGET_LINK_LIBRARIES(test_GainOperator diffusion utilities materials geometry)

ADD_EXECUTABLE(test_LossOperator        test_LossOperator.cc)
TARGET_LINK_LIBRARIES(test_LossOperator diffusion utilities materials geometry)

ADD_EXECUTABLE(test_OneGroupLossOperator        test_OneGroupLossOperator.cc)
TARGET_LINK_LIBRARIES(test_OneGroupLossOperator diffusion utilities materials geometry)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_GainOperator  	     test_GainOperator 0)
ADD_TEST( test_LossOperator  	     test_LossOperator 0)
ADD_TEST( test_OneGroupLossOperator  test_OneGroupLossOperator 0)
