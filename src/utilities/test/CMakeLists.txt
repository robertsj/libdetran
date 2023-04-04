#------------------------------------------------------------------------------#
# libdetran/src/utilities/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    utilities
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_Factory)
test(test_InputDB)
test(test_Iterators)
test(test_MathUtilities)
test(test_Random)
test(test_Testing)
test(test_TinyVector)