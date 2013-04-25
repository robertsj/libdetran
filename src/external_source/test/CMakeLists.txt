#------------------------------------------------------------------------------#
# src/transport/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_ConstantSource          test_ConstantSource.cc)
TARGET_LINK_LIBRARIES(test_ConstantSource   external_source)

ADD_EXECUTABLE(test_IsotropicSource         test_IsotropicSource.cc)
TARGET_LINK_LIBRARIES(test_IsotropicSource  external_source)

ADD_EXECUTABLE(test_DiscreteSource          test_DiscreteSource.cc)
TARGET_LINK_LIBRARIES(test_DiscreteSource   external_source)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_ConstantSource 		test_ConstantSource  0)
ADD_TEST(test_IsotropicSource      	test_IsotropicSource 0)
ADD_TEST(test_DiscreteSource       	test_DiscreteSource  0)

