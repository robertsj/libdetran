#------------------------------------------------------------------------------#
# libdetran/src/kinetics/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_KineticsMaterial            test_KineticsMaterial.cc)
TARGET_LINK_LIBRARIES(test_KineticsMaterial     kinetics)
ADD_EXECUTABLE(test_LinearMaterial              test_LinearMaterial.cc)
TARGET_LINK_LIBRARIES(test_LinearMaterial       kinetics)

ADD_EXECUTABLE(test_LinearExternalSource             test_LinearExternalSource.cc)
TARGET_LINK_LIBRARIES(test_LinearExternalSource      kinetics)
ADD_EXECUTABLE(test_PulsedExternalSource             test_PulsedExternalSource.cc)
TARGET_LINK_LIBRARIES(test_PulsedExternalSource      kinetics)

ADD_EXECUTABLE(test_MultiPhysics             		test_MultiPhysics.cc)
TARGET_LINK_LIBRARIES(test_MultiPhysics             kinetics)


#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_KineticsMaterial       test_KineticsMaterial 0)
ADD_TEST(test_LinearMaterial         test_LinearMaterial 0)

ADD_TEST(test_LinearExternalSource   test_LinearExternalSource 0)
ADD_TEST(test_PulsedExternalSource   test_PulsedExternalSource 0)

ADD_TEST(test_MultiPhysics           test_MultiPhysics 0)

