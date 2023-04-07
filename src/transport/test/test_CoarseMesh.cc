//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_CoarseMesh.cc
 *  @brief Test of CoarseMesh class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include <gtest/gtest.h>
#include "CoarseMesh.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"

#include "coarsemesh_fixture.hh"

using namespace detran;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;

TEST(CoarseMesh, Basic)
{
  //--------------------------------------------------------------------------//
  // SHARED DATA
  //--------------------------------------------------------------------------//

  // Coarse mesh dx's in all directions
  double delta[] = {3.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  //--------------------------------------------------------------------------//
  // 1D TESTS
  //--------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh1 = coarsemesh_1d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh1 = coarsemesh1->get_coarse_mesh();

  EXPECT_EQ(cmesh1->number_cells()  , 7);
  EXPECT_EQ(cmesh1->number_cells_x(), 7);
  EXPECT_EQ(cmesh1->number_cells_y(), 1);
  EXPECT_EQ(cmesh1->number_cells_z(), 1);
  for (int i = 0; i < cmesh1->number_cells(); i++)
  {
    EXPECT_NEAR(cmesh1->dx(i), delta[i], 1.0e-12);
  }
  EXPECT_NEAR(cmesh1->dy(0), 1.0, 1.0e-12);
  EXPECT_NEAR(cmesh1->dz(0), 1.0, 1.0e-12);

  const vec_int &mat1 = coarsemesh1->get_fine_mesh()->mesh_map("COARSEMESH");
  int refmat1[] = {0,0,0,1,1,2,2,3,3,4,4,5,5,6,6};
  for (int i = 0; i < mat1.size(); ++i)
  {
    EXPECT_EQ(refmat1[i], mat1[i]);
  }

  // Test for 1 coarse mesh.
  {
    vec_int fm(1, 5);
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int mt(1, 0);
    Mesh::SP_mesh mesh = std::make_shared<Mesh1D>(fm, cm, mt);
    CoarseMesh::SP_coarsemesh mesher(new CoarseMesh(mesh, 3));
    EXPECT_EQ(mesher->get_coarse_mesh()->number_cells(), 1);
  }

  //--------------------------------------------------------------------------//
  // 2D TESTS
  //--------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh2 = coarsemesh_2d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh2 = coarsemesh2->get_coarse_mesh();

  EXPECT_EQ(cmesh2->number_cells()  , 49);
  EXPECT_EQ(cmesh2->number_cells_x(), 7);
  EXPECT_EQ(cmesh2->number_cells_y(), 7);
  EXPECT_EQ(cmesh2->number_cells_z(), 1);
  for (int i = 0; i < cmesh2->number_cells_x(); i++)
  {
    EXPECT_NEAR(cmesh2->dx(i), delta[i], 1.0e-12);
    EXPECT_NEAR(cmesh2->dy(i), delta[i], 1.0e-12);
  }
  EXPECT_NEAR(cmesh2->dz(0), 1.0, 1.0e-12);

  //--------------------------------------------------------------------------//
  // 3D TESTS
  //--------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh3 = coarsemesh_3d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh3 = coarsemesh3->get_coarse_mesh();

  EXPECT_EQ(cmesh3->number_cells()  , 343);
  EXPECT_EQ(cmesh3->number_cells_x(), 7);
  EXPECT_EQ(cmesh3->number_cells_y(), 7);
  EXPECT_EQ(cmesh3->number_cells_z(), 7);
  const vec_int &mat_map3 = cmesh3->mesh_map("MATERIAL");
  for (int i = 0; i < cmesh3->number_cells_x(); i++)
  {
    EXPECT_NEAR(cmesh3->dx(i), delta[i], 1.0e-12);
    EXPECT_NEAR(cmesh3->dy(i), delta[i], 1.0e-12);
    EXPECT_NEAR(cmesh3->dz(i), delta[i], 1.0e-12);
    EXPECT_EQ(mat_map3[i], i);
  }
}

//----------------------------------------------------------------------------//
//              end of test_CoarseMesh.cc
//----------------------------------------------------------------------------//
