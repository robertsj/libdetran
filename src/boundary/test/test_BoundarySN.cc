//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  test_BoundarySN.cc
 *  @brief Test of BoundarySN class
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_BoundarySN)

#include "TestDriver.hh"
#include "boundary/BoundarySN.hh"
#include "angle/QuadratureFactory.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"

using namespace detran;
using namespace detran_geometry;
using namespace detran_angle;
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

int test_BoundarySN(int argc, char *argv[])
{

  // 1-D test
  {
    typedef BoundarySN<_1D> Boundary_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    inp->put<string>("bc_west",  "reflect");
    inp->put<string>("bc_left",  "vacuum");
    // Quadrature
    inp->put<int>("quad_number_polar_octant", 2);
    Boundary_T::SP_quadrature q = QuadratureFactory::build(inp, 1);
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh1D(fm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));
    // Tests
    TEST(b->is_reflective(0));
    TEST(!b->is_reflective(1));
    TEST(b->has_reflective());
    TEST(b->boundary_flux_size(0) == 4);
    TEST(b->boundary_flux_size(1) == 4);
  }

  // 2-D test
  {
    typedef BoundarySN<_2D> Boundary_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    inp->put<string>("bc_west",  "reflect");
    inp->put<string>("bc_left",  "vacuum");
    inp->put<string>("bc_south", "reflect");
    inp->put<string>("bc_north", "vacuum");
    // Quadrature
    inp->put<int>("quad_number_polar_octant", 2);
    inp->put<int>("quad_number_azimuth_octant", 2);
    Boundary_T::SP_quadrature q = QuadratureFactory::build(inp, 2);
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));
    // Tests
    TEST(b->is_reflective(0));
    TEST(!b->is_reflective(1));
    TEST(b->is_reflective(2));
    TEST(!b->is_reflective(3));
    TEST(b->has_reflective());
    TEST(b->boundary_flux_size(0) == 4*2*2*10);
    TEST(b->boundary_flux_size(1) == 4*2*2*10);
    TEST(b->boundary_flux_size(2) == 4*2*2*10);
    TEST(b->boundary_flux_size(3) == 4*2*2*10);
  }

  // 3-D test
  {
    typedef BoundarySN<_3D> Boundary_T;
    // Input
    Boundary_T::SP_input inp(new InputDB());
    inp->put<int>("number_groups", 2);
    inp->put<string>("bc_west",   "reflect");
    inp->put<string>("bc_left",   "vacuum");
    inp->put<string>("bc_south",  "reflect");
    inp->put<string>("bc_north",  "vacuum");
    inp->put<string>("bc_bottom", "reflect");
    inp->put<string>("bc_top",    "vacuum");
    // Quadrature
    inp->put<int>("quad_number_polar_octant",   2);
    inp->put<int>("quad_number_azimuth_octant", 2);
    Boundary_T::SP_quadrature q = QuadratureFactory::build(inp, 3);
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));
    // Tests
    TEST(b->is_reflective(0));
    TEST(!b->is_reflective(1));
    TEST(b->is_reflective(2));
    TEST(!b->is_reflective(3));
    TEST(b->is_reflective(4));
    TEST(!b->is_reflective(5));
    TEST(b->has_reflective());
    TEST(b->boundary_flux_size(0) == 8*2*2*10*10);
    TEST(b->boundary_flux_size(1) == 8*2*2*10*10);
    TEST(b->boundary_flux_size(2) == 8*2*2*10*10);
    TEST(b->boundary_flux_size(3) == 8*2*2*10*10);
    TEST(b->boundary_flux_size(4) == 8*2*2*10*10);
    TEST(b->boundary_flux_size(5) == 8*2*2*10*10);
  }

  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_BoundarySN.cc
//---------------------------------------------------------------------------//
