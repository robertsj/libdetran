#------------------------------------------------------------------------------#
# libdetran/src/geometry/test/CMakeLists.txt
#------------------------------------------------------------------------------#

macro(test)
  add_executable(${ARGV0} ${ARGV0}.cc)
  target_link_libraries(${ARGV0}
    geometry utilities angle
    GTest::gtest_main
  )
  gtest_discover_tests(${ARGV0})
endmacro()

test(test_CSG)
test(test_Geometry)
test(test_Mesh1D)
test(test_Mesh2D)
test(test_Mesh3D)
test(test_PinCell)
test(test_Point)
test(test_QuadraticSurface)
test(test_Ray)
test(test_Region)
test(test_RegionFactory)
test(test_Segment)
test(test_Track)
test(test_TrackDB)
test(test_Tracker)