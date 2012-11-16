#------------------------------------------------------------------------------#
# libdetran/src/kinetics/test/CMakeLists.txt
#------------------------------------------------------------------------------#

message( "++ Adding detran kinetics tests." )

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_KineticsMaterial            test_KineticsMaterial.cc)
TARGET_LINK_LIBRARIES(test_KineticsMaterial     kinetics)
ADD_EXECUTABLE(test_LinearMaterial              test_LinearMaterial.cc)
TARGET_LINK_LIBRARIES(test_LinearMaterial       kinetics)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_KineticsMaterial       test_KineticsMaterial 0)
ADD_TEST(test_LinearMaterial         test_LinearMaterial 0)
