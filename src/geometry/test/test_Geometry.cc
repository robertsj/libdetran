//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_Geometry.cc
 *  @brief Test of Geometry class
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "geometry/Geometry.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "callow/utils/Initialization.hh"
#include <cmath>
#define COUT(c) std::cout << c << std::endl;

using namespace detran_geometry;
using namespace detran_utilities;
using std::cout;
using std::endl;

TEST(Geometry, FromMesh2D)
{
  // make a 2-D mesh with dimensions and materials indices as follows
  //
  //   -------------------------   y=2.5
  //  |  3  |   2   |     1    |
  //  |-------------------------   y=2.0
  //  |     |       |          |
  //  |  0  |   1   |     2    |
  //  --------------------------   y=0.0
  // x=0.0  0.5     1.5        3.0

  std::vector<double> xfm {0.0, 0.5, 1.5, 3.0};
  std::vector<double> yfm {0.0, 2.0, 2.5};
  std::vector<int> mt { 0, 1, 2, 3, 2, 1};
  auto mesh = std::make_shared<Mesh2D>(xfm, yfm, mt);

  Geometry geo(mesh);

  EXPECT_EQ(geo.number_regions(), 6);

  // ensure proper bounding box
  EXPECT_EQ(geo.width_x(), 3.0);
  EXPECT_EQ(geo.width_y(), 2.5);
  EXPECT_EQ(geo.width_z(), 1.0);

  // ensure proper region ordering by getting region
  // for each mesh-cell center.
  EXPECT_EQ(geo.find(Point(0.25, 1.00)), 0);
  EXPECT_EQ(geo.find(Point(1.00, 1.00)), 1);
  EXPECT_EQ(geo.find(Point(2.25, 1.00)), 2);
  EXPECT_EQ(geo.find(Point(0.25, 2.25)), 3);
  EXPECT_EQ(geo.find(Point(1.00, 2.25)), 4);
  EXPECT_EQ(geo.find(Point(2.25, 2.25)), 5);

  // ensure proper material assignment
  EXPECT_EQ(geo.material_index(0), 0);
  EXPECT_EQ(geo.material_index(1), 1);
  EXPECT_EQ(geo.material_index(2), 2);
  EXPECT_EQ(geo.material_index(3), 3);
  EXPECT_EQ(geo.material_index(4), 2);
  EXPECT_EQ(geo.material_index(5), 1);
}

//----------------------------------------------------------------------------//
TEST(Geometry, FromMesh3D)
{
  // make a 3-D mesh with dimensions and materials indices as follows
  // so, in the top layer (z=1.0..2.0), the mt index is 10+ the mt of bottom.
  //       -----------------------  z=2.0
  //      /   13 /  12   /  11      /
  //     -----------------------  z=1.0
  //    /  3  /   2  /    1      /
  //   -------------------------  z=0.0,  y=2.5
  //  |  3  |   2   |     1    |
  //  |-------------------------   y=2.0
  //  |     |       |          |
  //  |  0  |   1   |     2    |
  //  --------------------------   y=0.0
  // x=0.0  0.5     1.5        3.0

  std::vector<double> xfm {0.0, 0.5, 1.5, 3.0};
  std::vector<double> yfm {0.0, 2.0, 2.5};
  std::vector<double> zfm {0.0, 1.0, 2.0};
  std::vector<int> mt { 0, 1, 2, 3, 2, 1,
                       10,11,12,13,12,11};
  auto mesh = std::make_shared<Mesh3D>(xfm, yfm, zfm, mt);

  Geometry geo(mesh);

  EXPECT_EQ(geo.number_regions(), 12);

  // ensure proper bounding box
  EXPECT_EQ(geo.width_x(), 3.0);
  EXPECT_EQ(geo.width_y(), 2.5);
  EXPECT_EQ(geo.width_z(), 2.0);

  // ensure proper region ordering by getting region
  // for each mesh-cell center.
  EXPECT_EQ(geo.find(Point(0.25, 1.00, 0.5)), 0);
  EXPECT_EQ(geo.find(Point(1.00, 1.00, 0.5)), 1);
  EXPECT_EQ(geo.find(Point(2.25, 1.00, 0.5)), 2);
  EXPECT_EQ(geo.find(Point(0.25, 2.25, 0.5)), 3);
  EXPECT_EQ(geo.find(Point(1.00, 2.25, 0.5)), 4);
  EXPECT_EQ(geo.find(Point(2.25, 2.25, 0.5)), 5);
  EXPECT_EQ(geo.find(Point(0.25, 1.00, 1.5)), 6);
  EXPECT_EQ(geo.find(Point(1.00, 1.00, 1.5)), 7);
  EXPECT_EQ(geo.find(Point(2.25, 1.00, 1.5)), 8);
  EXPECT_EQ(geo.find(Point(0.25, 2.25, 1.5)), 9);
  EXPECT_EQ(geo.find(Point(1.00, 2.25, 1.5)), 10);
  EXPECT_EQ(geo.find(Point(2.25, 2.25, 1.5)), 11);

  // ensure proper material assignment
  EXPECT_EQ(geo.material_index(0), 0);
  EXPECT_EQ(geo.material_index(1), 1);
  EXPECT_EQ(geo.material_index(2), 2);
  EXPECT_EQ(geo.material_index(3), 3);
  EXPECT_EQ(geo.material_index(4), 2);
  EXPECT_EQ(geo.material_index(5), 1);
  EXPECT_EQ(geo.material_index(6), 10);
  EXPECT_EQ(geo.material_index(7), 11);
  EXPECT_EQ(geo.material_index(8), 12);
  EXPECT_EQ(geo.material_index(9), 13);
  EXPECT_EQ(geo.material_index(10), 12);
  EXPECT_EQ(geo.material_index(11), 11);
}

//----------------------------------------------------------------------------//
//              end of test_Geometry.cc
//----------------------------------------------------------------------------//
