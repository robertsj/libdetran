#------------------------------------------------------------------------------#
# src/transport/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    external_source
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_ConstantSource)
test(test_DiscreteSource)
test(test_IsotropicSource)
