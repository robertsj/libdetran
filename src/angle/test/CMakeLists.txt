#------------------------------------------------------------------------------#
# libdetran/src/angle/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    angle
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_BaseQuadrature)
#test(test_Collocated)
test(test_MomentIndexer)
test(test_MomentToDiscrete)
test(test_PolarQuadrature)
test(test_ProductQuadratures)
#test(test_QuadrupleRange)
test(test_SphericalHarmonics)
#test(test_TabuchiYamamoto)
#test(test_Uniform)