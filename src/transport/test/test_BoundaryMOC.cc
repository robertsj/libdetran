//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_BoundaryMOC.cc
 * \author Jeremy Roberts
 * \date   Jul 17, 2012
 * \brief  Test of BoundaryMOC.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_BoundaryMOC)        \
        FUNC(test_BoundaryMOC_2)

// Detran headers
#include "TestDriver.hh"
#include "BoundaryMOC.hh"
#include "Mesh2D.hh"
#include "Collocated.hh"
#include "Uniform.hh"
#include "Tracker.hh"

// Setup
#include "geometry/test/mesh_fixture.hh"
#include "angle/test/quadrature_fixture.hh"

using namespace detran;
using namespace detran_test;
using namespace std;

int main(int argc, char *argv[])
{
  RUN(argc, argv);
}

//----------------------------------------------//
// TEST DEFINITIONS
//----------------------------------------------//

int test_BoundaryMOC()
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh0(new Mesh2D(fm, fm, cm, cm, mt));

  // Create quadrature.  1 azimuth, 3 space, 1 polar.
  QuadratureMOC::SP_quadrature quad(new Uniform(2, 2, 5, 3, "TY"));

  // Create tracker
  Tracker tracker(mesh0, quad);
  Tracker::SP_trackdb tracks = tracker.trackdb();

  // Get the tracked mesh.
  Mesh::SP_mesh mesh = tracker.meshmoc();

  // Empty input
  InputDB::SP_input input(new InputDB());
  input->put<int>("number_groups", 1);

  // Create the boundary
  BoundaryMOC<_2D> boundary(input, mesh, quad);

  // Now test things.  Four octants, one angle each, three tracks
  // each.  The tracks are paired as follows:
  //
  // (0,0,0) -> (3,0,0)
  // (0,0,1) -> (1,0,0)
  // (0,0,2) -> (1,0,1)
  //
  // (1,0,0) -> (2,0,0)
  // (1,0,1) -> (0,0,0)
  // (1,0,2) -> (0,0,1)
  //
  // (2,0,0) -> (3,0,1)
  // (2,0,1) -> (3,0,2)
  // (2,0,2) -> (1,0,2)
  //
  // (3,0,0) -> (2,0,1)
  // (3,0,1) -> (2,0,2)
  // (3,0,2) -> (0,0,2)
  //

  detran::u_int o;
  detran::u_int a;
  detran::u_int t;

  if (1==0)
  {

  // Octant 0
  boundary.feed_into(0, 0, 0, o, a, t);
  TEST(o == 3);
  TEST(t == 0);
  boundary.feed_into(0, 0, 1, o, a, t);
  TEST(o == 1);
  TEST(t == 0);
  boundary.feed_into(0, 0, 2, o, a, t);
  TEST(o == 1);
  TEST(t == 1);

  // Octant 1
  boundary.feed_into(1, 0, 0, o, a, t);
  TEST(o == 2);
  TEST(t == 0);
  boundary.feed_into(1, 0, 1, o, a, t);
  TEST(o == 0);
  TEST(t == 0);
  boundary.feed_into(1, 0, 2, o, a, t);
  TEST(o == 0);
  TEST(t == 1);

  // Octant 2
  boundary.feed_into(2, 0, 0, o, a, t);
  TEST(o == 3);
  TEST(t == 1);
  boundary.feed_into(2, 0, 1, o, a, t);
  TEST(o == 3);
  TEST(t == 2);
  boundary.feed_into(2, 0, 2, o, a, t);
  TEST(o == 1);
  TEST(t == 2);

  // Octant 3
  boundary.feed_into(3, 0, 0, o, a, t);
  TEST(o == 2);
  TEST(t == 1);
  boundary.feed_into(3, 0, 1, o, a, t);
  TEST(o == 2);
  TEST(t == 2);
  boundary.feed_into(3, 0, 2, o, a, t);
  TEST(o == 0);
  TEST(t == 2);
  }


  // Test side indices
  for (int s = 0; s < 4; s++)
  {
    cout << " side = " << s << endl;
    vec2_int left = boundary.side_indices(s);
    for (int i = 0; i < left.size(); i++)
    {
      cout << left[i][0] << " " << left[i][1] << " " << left[i][2] << " | ";

      detran::u_int o = left[i][0];
      detran::u_int a = left[i][1];
      detran::u_int t = left[i][2];
      detran::u_int oo;
      detran::u_int aa;
      detran::u_int tt;
      boundary.feed_into(o,a,t,oo,aa,tt);
//      cout << " INTO " << endl;
//      cout << "      o = " << oo << endl
//           << "        a = " << aa << endl
//           << "          t = " << tt << endl;
      cout << oo << " " << aa << " " << tt << " | ";

      boundary.feed_from(o,a,t,oo,aa,tt);
//      cout << " FROM " << endl;
//      cout << "      o = " << oo << endl
//           << "        a = " << aa << endl
//           << "          t = " << tt << endl;
      cout << oo << " " << aa << " " << tt << endl;

    }
  }


  return 0;
}

int test_BoundaryMOC_2()
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mt(1, 0);
  Mesh::SP_mesh mesh0(new Mesh2D(fm, fm, cm, cm, mt));

  // Create quadrature.  1 azimuth, 3 space, 1 polar.
  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));

  // Create tracker
  Tracker tracker(mesh0, quad);
  Tracker::SP_trackdb tracks = tracker.trackdb();

  // Get the tracked mesh.
  Mesh::SP_mesh mesh = tracker.meshmoc();

  // Empty input
  InputDB::SP_input input(new InputDB());
  input->put<int>("number_groups", 1);

  // Create the boundary
  BoundaryMOC<_2D> boundary(input, mesh, quad);

  // Now test things.  Four octants, one angle each, three tracks
  // each.  The tracks are paired as follows:
  //
  // (0,0, 0) -> (3,0,0)
  // (0,0, 1) -> (1,0,0)
  // (0,0, 2) -> (1,0,1)
  // (0,0, 3) -> (3,0,0)
  // (0,0, 4) -> (1,0,0)
  // (0,0, 5) -> (1,0,1)
  // (0,0, 6) -> (3,0,0)
  // (0,0, 7) -> (1,0,0)
  // (0,0, 8) -> (1,0,1)
  // (0,0, 9) -> (3,0,0)
  // (0,0,10) -> (1,0,0)
  // (0,0,11) -> (1,0,1)
  // (0,0,12) -> (3,0,0)
  // (0,0,13) -> (1,0,0)
  //
  // (2,0, 0) -> (3,0,0)
  // (2,0, 1) -> (1,0,0)
  // (2,0, 2) -> (1,0,1)
  // (2,0, 3) -> (3,0,0)
  // (2,0, 4) -> (1,0,0)
  // (2,0, 5) -> (1,0,1)
  // (2,0, 6) -> (3,0,0)
  // (2,0, 7) -> (1,0,0)
  // (2,0, 8) -> (1,0,1)
  // (2,0, 9) -> (3,0,0)
  // (2,0,10) -> (1,0,0)
  // (2,0,11) -> (1,0,1)
  // (2,0,12) -> (3,0,0)
  // (2,0,13) -> (1,0,0)
  //
  detran::u_int o;
  detran::u_int a;
  detran::u_int t;

  // Octant 0
  boundary.feed_into(0, 0, 0, o, a, t);
  TEST(o == 3);
  TEST(t == 0);
  boundary.feed_into(0, 0, 1, o, a, t);
  TEST(o == 1);
  TEST(t == 0);
  boundary.feed_into(0, 0, 2, o, a, t);
  TEST(o == 1);
  TEST(t == 1);

  // Octant 1
  boundary.feed_into(1, 0, 0, o, a, t);
  TEST(o == 2);
  TEST(t == 0);
  boundary.feed_into(1, 0, 1, o, a, t);
  TEST(o == 0);
  TEST(t == 0);
  boundary.feed_into(1, 0, 2, o, a, t);
  TEST(o == 0);
  TEST(t == 1);

  // Octant 2
  boundary.feed_into(2, 0, 0, o, a, t);
  TEST(o == 3);
  TEST(t == 1);
  boundary.feed_into(2, 0, 1, o, a, t);
  TEST(o == 3);
  TEST(t == 2);
  boundary.feed_into(2, 0, 2, o, a, t);
  TEST(o == 1);
  TEST(t == 2);

  // Octant 3
  boundary.feed_into(3, 0, 0, o, a, t);
  TEST(o == 2);
  TEST(t == 1);
  boundary.feed_into(3, 0, 1, o, a, t);
  TEST(o == 2);
  TEST(t == 2);
  boundary.feed_into(3, 0, 2, o, a, t);
  TEST(o == 0);
  TEST(t == 2);

  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_State.cc
//---------------------------------------------------------------------------//
