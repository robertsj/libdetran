#------------------------------------------------------------------------------#
# libdetran/src/angle/test/CMakeLists.txt
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# Create test executables
#------------------------------------------------------------------------------#

ADD_EXECUTABLE(test_MomentToDiscrete        test_MomentToDiscrete.cc)
TARGET_LINK_LIBRARIES(test_MomentToDiscrete angle)

ADD_EXECUTABLE(test_MomentIndexer           test_MomentIndexer.cc)
TARGET_LINK_LIBRARIES(test_MomentIndexer    angle)
ADD_TEST(test_MomentIndexer_1D              test_MomentIndexer      0)
ADD_TEST(test_MomentIndexer_1D              test_MomentIndexer      1)
ADD_TEST(test_MomentIndexer_1D              test_MomentIndexer      2)

ADD_EXECUTABLE(test_SphericalHarmonics        test_SphericalHarmonics.cc)
TARGET_LINK_LIBRARIES(test_SphericalHarmonics angle)

ADD_EXECUTABLE(test_ProductQuadratures        test_ProductQuadratures.cc)
TARGET_LINK_LIBRARIES(test_ProductQuadratures angle)

ADD_EXECUTABLE(test_BaseQuadrature          test_BaseQuadrature.cc)
TARGET_LINK_LIBRARIES(test_BaseQuadrature   angle)
ADD_TEST(test_BaseQuadrature_GL             test_BaseQuadrature      0)

ADD_EXECUTABLE(test_PolarQuadrature         test_PolarQuadrature.cc)
TARGET_LINK_LIBRARIES(test_PolarQuadrature  angle)
ADD_TEST(test_PolarQuadrature_basic         test_PolarQuadrature      0)

#------------------------------------------------------------------------------#
# Define tests
#------------------------------------------------------------------------------#

ADD_TEST(test_MomentToDiscrete      test_MomentToDiscrete   0)
ADD_TEST(test_SphericalHarmonics    test_SphericalHarmonics 0)
