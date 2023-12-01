#------------------------------------------------------------------------------#
# src/boundary/test/CMakeLists.txt
#------------------------------------------------------------------------------#
macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    boundary
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_BoundaryDiffusion)
test(test_BoundaryMOC)
test(test_BoundarySN)