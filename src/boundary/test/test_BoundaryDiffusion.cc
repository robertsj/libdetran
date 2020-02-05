//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_BoundaryDiffusion.cc
 *  @brief Test of BoundaryDiffusion class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                       \
        FUNC(test_BoundaryDiffusion_1D) \
        FUNC(test_BoundaryDiffusion_2D) \
        FUNC(test_BoundaryDiffusion_3D)

#include "utilities/TestDriver.hh"
#include "boundary/BoundaryDiffusion.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"

using namespace detran;
using namespace detran_geometry;
using namespace detran_utilities;
using namespace detran_test;
using std::cout;
using std::endl;
using std::string;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------------------------------------//

int test_BoundaryDiffusion_1D(int argc, char *argv[])
{
  // 1-D test
  {
    typedef BoundaryDiffusion<_1D> Boundary_T;
    typedef BoundaryValue<_1D>     BV_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh1D(fm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh));
    // Fill the boundary
    Boundary_T &b_ref = *b;
    cout << " bref = " << b_ref(Mesh::WEST, 0, Boundary_T::IN) << endl;
    (*b)(Mesh::WEST, 0, Boundary_T::IN) = 1.0;
    (*b)(Mesh::WEST, 1, Boundary_T::IN) = 2.0;
    (*b)(Mesh::EAST, 0, Boundary_T::IN) = 3.0;
    (*b)(Mesh::EAST, 1, Boundary_T::IN) = 4.0;
    (*b)(Mesh::WEST, 0, Boundary_T::OUT) = 1.1;
    (*b)(Mesh::WEST, 1, Boundary_T::OUT) = 2.2;
    (*b)(Mesh::EAST, 0, Boundary_T::OUT) = 3.3;
    (*b)(Mesh::EAST, 1, Boundary_T::OUT) = 4.4;
    // Test
    TEST(soft_equiv( BV_T::value((*b)(Mesh::WEST, 0, Boundary_T::IN)), 1.0 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::WEST, 1, Boundary_T::IN)), 2.0 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::EAST, 0, Boundary_T::IN)), 3.0 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::EAST, 1, Boundary_T::IN)), 4.0 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::WEST, 0, Boundary_T::OUT)), 1.1 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::WEST, 1, Boundary_T::OUT)), 2.2 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::EAST, 0, Boundary_T::OUT)), 3.3 ));
    TEST(soft_equiv( BV_T::value((*b)(Mesh::EAST, 1, Boundary_T::OUT)), 4.4 ));
  }
  return 0;
}

int test_BoundaryDiffusion_2D(int argc, char *argv[])
{
  // 2-D test
  {
    typedef BoundaryDiffusion<_2D> Boundary_T;
    typedef BoundaryValue<_2D>     BV_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    inp->put<string>("bc_west",  "reflect");
    inp->put<string>("bc_left",  "vacuum");
    inp->put<string>("bc_south", "reflect");
    inp->put<string>("bc_north", "vacuum");
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
    // Boundary
    Boundary_T b(inp, mesh);
    for (int s = 0; s < 4; ++s)
    {
      for (int i = 0; i < 10; ++i)
      {
        BV_T::value(b(s, 0, Boundary_T::IN),  i) = 100 * s + i;
        BV_T::value(b(s, 1, Boundary_T::IN),  i) = 100 * s + i + 1000;
        BV_T::value(b(s, 0, Boundary_T::OUT), i) = 100 * s + 0.1*i;
        BV_T::value(b(s, 1, Boundary_T::OUT), i) = 100 * s + 0.1*i + 1000;
      }
    }
    // Test
    TEST(soft_equiv( BV_T::value(b(Mesh::WEST,  0, Boundary_T::IN),  5), 5.0));
    TEST(soft_equiv( BV_T::value(b(Mesh::NORTH, 1, Boundary_T::OUT), 2), 1300.2));
  }
  return 0;
}

int test_BoundaryDiffusion_3D(int argc, char *argv[])
{
  // 3-D test
  {
    typedef BoundaryDiffusion<_3D> Boundary_T;
    typedef BoundaryValue<_3D>     BV_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    inp->put<string>("bc_west",   "reflect");
    inp->put<string>("bc_left",   "vacuum");
    inp->put<string>("bc_south",  "reflect");
    inp->put<string>("bc_north",  "vacuum");
    inp->put<string>("bc_bottom", "reflect");
    inp->put<string>("bc_top",    "vacuum");
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));
    // Boundary
    Boundary_T b(inp, mesh);
    for (int s = 0; s < 6; ++s)
    {
      for (int i = 0; i < 10; ++i)
      {
        for (int j = 0; j < 10; ++j)
        {
          BV_T::value(b(s, 0, Boundary_T::IN),  i, j) = 100 * s + i + 10*j;
          BV_T::value(b(s, 1, Boundary_T::IN),  i, j) = 100 * s + i + 10*j + 1000;
          BV_T::value(b(s, 0, Boundary_T::OUT), i, j) = 100 * s + 0.1*i + 10*j;
          BV_T::value(b(s, 1, Boundary_T::OUT), i, j) = 100 * s + 0.1*i + 10*j + 1000;
        }
      }
    }
    // Test
    TEST(soft_equiv( BV_T::value(b(Mesh::WEST,  0, Boundary_T::IN),  5, 5),   55.0));
    TEST(soft_equiv( BV_T::value(b(Mesh::NORTH, 1, Boundary_T::OUT), 2, 2), 1320.2));
  }
  return 0;
}

//----------------------------------------------------------------------------//
//              end of test_BoundaryDiffusion.cc
//----------------------------------------------------------------------------//
