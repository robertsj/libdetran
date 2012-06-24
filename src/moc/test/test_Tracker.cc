//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test_Tracker.cc
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 * \brief  Test of Tracker class
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// LIST OF TEST FUNCTIONS
#define TEST_LIST                     \
        FUNC(test_Tracker_2x2)        \
        FUNC(test_Tracker_3x3)

// Detran headers
#include "TestDriver.hh"
#include "Tracker.hh"
//
#include "Mesh2D.hh"
#include "Uniform.hh"

// Setup
/* ... */

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

int test_Tracker_2x2()
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 2);
  vec_int mat(1, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mat));

  // Create quadrature
  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
  quad->display_tracks();

  // Create tracker
  Tracker tracker(mesh, quad);
  Tracker::SP_trackdb tracks = tracker.trackdb();

  // Verify tracker.
  double length = 0.559016994374947;
  // Number of segments for each track within angle
  int ns[] = {1,2,1};
  // Region map for all segments
  int region[] = {2,0,3,1,3,1,2,0};
  // Region counter
  int r = 0;
  for (int a = 0; a < 2; a++)
  {
    TEST(tracks->number_tracks_angle(a) == 3);
    for (int t = 0; t < 3; t++)
    {
      Tracker::SP_track track = tracks->track(a, t);
      TEST(track->number_segments() == ns[t]);
      for (int s = 0; s < ns[t]; s++)
      {
        TEST(soft_equiv(track->segment(s).length(), length));
        TEST(track->segment(s).region() == region[r++]);
      }
    }
  }
  return 0;
}

int test_Tracker_3x3()
{

  // Create mesh
  vec_dbl cm(2, 0.0);
  cm[1] = 1.0;
  vec_int fm(1, 3);
  vec_int mat(1, 0);
  Mesh::SP_mesh mesh(new Mesh2D(fm, fm, cm, cm, mat));
  mesh->display();
  // Create quadrature
  QuadratureMOC::SP_quadrature quad(new Uniform(2, 1, 3, 1, "TY"));
  quad->display_tracks();

  // Create tracker
  Tracker tracker(mesh, quad);
  Tracker::SP_trackdb tracks = tracker.trackdb();
  tracks->display();

  // Number of segments for each track within angle
  int ns[] = {2,5,2};
  // Region map for all segments
  int region[] = {6,7, 0,3,4,5,8, 1,2,
                  8,7, 2,5,4,3,6, 1,0};
  double len0 = 0.186338998124982;
  double len1 = 0.372677996249965;
  double length[] = {len1,len0,  len0,len0,len1,len0,len0,  len0,len1};

  // Region counter
  int r = 0;
  for (int a = 0; a < 2; a++)
  {
    TEST(tracks->number_tracks_angle(a) == 3);
    // Length counter
    int l = 0;
    for (int t = 0; t < 3; t++)
    {
      Tracker::SP_track track = tracks->track(a, t);
      TEST(track->number_segments() == ns[t]);
      for (int s = 0; s < ns[t]; s++)
      {
        TEST(soft_equiv(track->segment(s).length(), length[l++]));
        TEST(track->segment(s).region() == region[r++]);
      }
    }
  }
  return 0;
}

//---------------------------------------------------------------------------//
//              end of test_Tracker.cc
//---------------------------------------------------------------------------//
