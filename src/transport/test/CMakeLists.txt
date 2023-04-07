#------------------------------------------------------------------------------#
# src/transport/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    transport
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_CoarseMesh)
test(test_CurrentTally)
test(test_Equation_DD_1D)
test(test_Equation_DD_2D)
test(test_Equation_SC_1D)
#test(test_Equation_SC_MOC)
test(test_Homogenization)
test(test_ScatterSource)
test(test_State)
test(test_Sweeper2D)
test(test_Sweeper3D)