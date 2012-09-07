//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_BoundarySN.cc
 * \author Jeremy Roberts
 * \date   Apr 1, 2012
 * \brief  Test of BoundarySN class
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                    \
        FUNC(test_BoundarySN)

// Detran headers
#include "TestDriver.hh"
#include "boundary/BoundarySN.hh"
#include "angle/GaussLegendre.hh"
#include "angle/QuadrupleRange.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"

// Setup
/* ... */

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

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

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
    Boundary_T::SP_quadrature q(new GaussLegendre(2));
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh1D(fm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));

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
    Boundary_T::SP_quadrature q(new QuadrupleRange(2, 2));
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));

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
    Boundary_T::SP_quadrature q(new QuadrupleRange(2, 3));
    // Mesh
    vec_dbl cm(2, 0.0); cm[1] = 1.0;
    vec_int fm(1, 10);
    vec_int mt(1, 0);
    Boundary_T::SP_mesh mesh(new Mesh3D(fm, fm, fm, cm, cm, cm, mt));
    // Boundary
    Boundary_T::SP_boundary b(new Boundary_T(inp, mesh, q));

  }


  return 0;
}


//---------------------------------------------------------------------------//
//              end of test_BoundarySN.cc
//---------------------------------------------------------------------------//
