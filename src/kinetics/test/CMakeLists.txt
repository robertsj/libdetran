#------------------------------------------------------------------------------#
# libdetran/src/kinetics/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    kinetics
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_KineticsMaterial)
test(test_LinearExternalSource)
test(test_LinearMaterial)
test(test_MultiPhysics)
test(test_PulsedExternalSource)
test(test_TimeIndependentMaterial)