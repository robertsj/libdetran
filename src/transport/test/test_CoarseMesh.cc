//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_CoarseMesh.cc
 *  @brief Test of CoarseMesh class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_CoarseMesh)

#include "utilities/TestDriver.hh"
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

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_CoarseMesh(int argc, char *argv[])
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

  TEST(cmesh1->number_cells()   == 7);
  TEST(cmesh1->number_cells_x() == 7);
  TEST(cmesh1->number_cells_y() == 1);
  TEST(cmesh1->number_cells_z() == 1);
  for (int i = 0; i < cmesh1->number_cells(); i++)
  {
    TEST(soft_equiv(cmesh1->dx(i), delta[i]));
  }
  TEST(soft_equiv(cmesh1->dy(0), 1.0));
  TEST(soft_equiv(cmesh1->dz(0), 1.0));

  const vec_int &mat1 = coarsemesh1->get_fine_mesh()->mesh_map("COARSEMESH");
  int refmat1[] = {0,0,0,1,1,2,2,3,3,4,4,5,5,6,6};
  for (int i = 0; i < mat1.size(); ++i)
  {
    TEST(refmat1[i] == mat1[i]);
  }

  // Test for 1 coarse mesh.
  {
    vec_int fm(1, 5);
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int mt(1, 0);
    Mesh::SP_mesh mesh = Mesh1D::Create(fm, cm, mt);
    CoarseMesh::SP_coarsemesh mesher(new CoarseMesh(mesh, 3));
    TEST(mesher->get_coarse_mesh()->number_cells() == 1);
  }

  //--------------------------------------------------------------------------//
  // 2D TESTS
  //--------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh2 = coarsemesh_2d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh2 = coarsemesh2->get_coarse_mesh();

  TEST(cmesh2->number_cells()   == 49);
  TEST(cmesh2->number_cells_x() == 7);
  TEST(cmesh2->number_cells_y() == 7);
  TEST(cmesh2->number_cells_z() == 1);
  for (int i = 0; i < cmesh2->number_cells_x(); i++)
  {
    TEST(soft_equiv(cmesh2->dx(i), delta[i]));
    TEST(soft_equiv(cmesh2->dy(i), delta[i]));
  }
  TEST(soft_equiv(cmesh2->dz(0), 1.0));

  //--------------------------------------------------------------------------//
  // 3D TESTS
  //--------------------------------------------------------------------------//

  // Get the coarsemesh from the fixture.
  CoarseMesh::SP_coarsemesh coarsemesh3 = coarsemesh_3d();

  // Get the underlying coarse mesh.
  Mesh::SP_mesh cmesh3 = coarsemesh3->get_coarse_mesh();

  TEST(cmesh3->number_cells()   == 343);
  TEST(cmesh3->number_cells_x() == 7);
  TEST(cmesh3->number_cells_y() == 7);
  TEST(cmesh3->number_cells_z() == 7);
  const vec_int &mat_map3 = cmesh3->mesh_map("MATERIAL");
  for (int i = 0; i < cmesh3->number_cells_x(); i++)
  {
    TEST(soft_equiv(cmesh3->dx(i), delta[i]));
    TEST(soft_equiv(cmesh3->dy(i), delta[i]));
    TEST(soft_equiv(cmesh3->dz(i), delta[i]));
    TEST(mat_map3[i] == i);
  }

  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_CoarseMesh.cc
//----------------------------------------------------------------------------//
