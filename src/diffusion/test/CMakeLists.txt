#------------------------------------------------------------------------------#
# /src/diffusion/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran diffusion tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_OneGroupLossOperator        test_OneGroupLossOperator.cc)
TARGET_LINK_LIBRARIES(test_OneGroupLossOperator diffusion utilities materials geometry)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST( test_OneGroupLossOperator  test_OneGroupLossOperator 0)
