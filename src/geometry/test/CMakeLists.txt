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

test(test_Point)
test(test_Tracker)

ADD_EXECUTABLE(test_Mesh1D                  test_Mesh1D.cc)
TARGET_LINK_LIBRARIES(test_Mesh1D           geometry utilities angle)
ADD_TEST(test_Mesh1D                        test_Mesh1D 0)

ADD_EXECUTABLE(test_Mesh2D                  test_Mesh2D.cc)
TARGET_LINK_LIBRARIES(test_Mesh2D           geometry utilities angle)
ADD_TEST(test_Mesh2D                        test_Mesh2D 0)

ADD_EXECUTABLE(test_Mesh3D                  test_Mesh3D.cc)
TARGET_LINK_LIBRARIES(test_Mesh3D           geometry utilities angle)
ADD_TEST(test_Mesh3D                        test_Mesh3D 0)

ADD_EXECUTABLE(test_PinCell                 test_PinCell.cc)
TARGET_LINK_LIBRARIES(test_PinCell          geometry utilities angle)
ADD_TEST(test_PinCell                       test_PinCell 0)

ADD_EXECUTABLE(test_Segment                 test_Segment.cc)
TARGET_LINK_LIBRARIES(test_Segment          geometry utilities angle)
ADD_TEST(test_Segment                       test_Segment 0)

ADD_EXECUTABLE(test_Track                   test_Track.cc)
TARGET_LINK_LIBRARIES(test_Track            geometry utilities angle)
ADD_TEST(test_Track                         test_Track 0)

ADD_EXECUTABLE(test_TrackDB                 test_TrackDB.cc)
TARGET_LINK_LIBRARIES(test_TrackDB          geometry utilities angle)
ADD_TEST(test_TrackDB_2D                    test_TrackDB 0)
ADD_TEST(test_TrackDB_3D                    test_TrackDB 1)

ADD_EXECUTABLE(test_QuadraticSurface        test_QuadraticSurface.cc)
TARGET_LINK_LIBRARIES(test_QuadraticSurface geometry utilities angle)
ADD_TEST(test_QuadraticSurface_plane        test_QuadraticSurface 0)
ADD_TEST(test_QuadraticSurface_circle       test_QuadraticSurface 1)

ADD_EXECUTABLE(test_Ray                     test_Ray.cc)
TARGET_LINK_LIBRARIES(test_Ray              geometry utilities angle)
ADD_TEST(test_Ray                           test_Ray 0)

ADD_EXECUTABLE(test_CSG                     test_CSG.cc)
TARGET_LINK_LIBRARIES(test_CSG              geometry utilities angle)
ADD_TEST(test_CSG_Intersection              test_CSG 0)
ADD_TEST(test_CSG_Union                     test_CSG 1)
ADD_TEST(test_CSG_Difference                test_CSG 2)
ADD_TEST(test_CSG_Translation               test_CSG 3)

ADD_EXECUTABLE(test_RegionFactory           test_RegionFactory.cc)
TARGET_LINK_LIBRARIES(test_RegionFactory    geometry utilities angle)
ADD_TEST(test_RegionFactory_PinCell         test_RegionFactory 0)

ADD_EXECUTABLE(test_Region                  test_Region.cc)
TARGET_LINK_LIBRARIES(test_Region           geometry utilities angle)
ADD_TEST(test_Region                        test_Region 0)

ADD_EXECUTABLE(test_Geometry                test_Geometry.cc)
TARGET_LINK_LIBRARIES(test_Geometry         geometry utilities angle)
ADD_TEST(test_Geometry                      test_Geometry 0)


