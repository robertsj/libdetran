macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    material utilities
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_Material)
