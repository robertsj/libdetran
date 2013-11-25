//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  csg_fixture.hh
 *  @brief Several CSG fixtures for testing
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_CSG_FIXTURE_HH_
#define detran_geometry_CSG_FIXTURE_HH_

#include "geometry/Geometry.hh"
#include "geometry/QuadraticSurfaceFactory.hh"
#include "geometry/RegionFactory.hh"
#include "utilities/MathUtilities.hh"

namespace detran_geometry
{

/// A 2-D pincell with two regions
static Geometry::SP_geometry test_2D_pincell_simple()
{
  typedef Surface::SP_surface     SP_surface;
  typedef Region::SP_region       SP_region;
  typedef QuadraticSurfaceFactory QSF;

  double P  = 1.26;    // pitch
  double HP = 0.5 * P; // half pitch
  double R  = 0.54;    // pin radius

  Geometry::SP_geometry geo = Geometry::Create(P, P, 0.0);

  // Surfaces
  SP_surface W = QSF::CreatePlaneX(0.0);
  SP_surface E = QSF::CreatePlaneX(P);
  SP_surface S = QSF::CreatePlaneY(0.0);
  SP_surface N = QSF::CreatePlaneY(P);
  SP_surface C = QSF::CreateCylinderZ(HP, HP, R);

  // Pin
  SP_region pin = Region::Create(0, Point(0,0,0)-1.0e-8, Point(P,P,0)+1.0e-8);
  pin->append(C, false);
  geo->add_region(pin);

  // Moderator
  SP_region mod = Region::Create(1, Point(0,0,0)-1.0e-8, Point(P,P,0)+1.0e-8);
  mod->append(C, true);
  mod->append(W, true);
  mod->append(E, false);
  mod->append(S, true);
  mod->append(N, false);
  geo->add_region(mod);

  return geo;
}

/// A 2-D pincell with slightly more complex division
static Geometry::SP_geometry test_2D_pincell_complex()
{
  typedef Surface::SP_surface     SP_surface;
  typedef Region::SP_region       SP_region;
  typedef QuadraticSurfaceFactory QSF;

  double P  = 1.26;    // pitch
  double HP = 0.5 * P; // half pitch
  double R  = 0.54;    // pin radius

  Geometry::SP_geometry geo = Geometry::Create(P, P, 0.0);

  // Surfaces
  SP_surface surfaces[] = {QSF::CreatePlaneX(0.0),                 // 1  W
                           QSF::CreatePlaneX(P),                   // 2  E
                           QSF::CreatePlaneY(0.0),                 // 3  S
                           QSF::CreatePlaneY(P),                   // 4  N
                           QSF::CreatePlaneX(HP),                  // 5  V
                           QSF::CreatePlaneY(HP),                  // 6  H
                           QSF::CreateCylinderZ(HP, HP, R),        // 7  O
                           QSF::CreateCylinderZ(HP, HP, 0.5 * R)}; // 8  I

  // Region surface lists
  int slist[12][5] = {{-8,  6,  5,  9,  9}, // inner, I
                      {-8,  6, -5,  9,  9}, // inner, II
                      {-8, -6, -5,  9,  9}, // inner, III
                      {-8, -6,  5,  9,  9}, // inner, IV
                      {-7,  8,  5,  6,  9}, // outer, I
                      {-7,  8,  5, -6,  9}, // outer, II
                      {-7,  8, -5, -6,  9}, // outer, III
                      {-7,  8, -5,  6,  9}, // outer, IV
                      { 7,  6, -4, -2,  5}, // mod, I
                      { 7,  6, -4,  1, -5}, // mod, II
                      { 7, -6,  3,  1, -5}, // mod, III
                      { 7, -6,  3, -2,  5}};// mod, IV

  // Region material map
  int mat[12] = {0,0,0,0,1,1,1,1,2,2,2,2};

  // Build the geometry
  for (int r = 0; r < 12; ++r)
  {
    SP_region tmp =
      Region::Create(mat[r], Point(0,0,0)-1.0e-8, Point(P,P,0)+1.0e-8);
    for (int s = 0; s < 5; ++s)
    {
      if (slist[r][s] == 9) break;
      tmp->append(surfaces[std::abs(slist[r][s])-1], slist[r][s] > 0);
    }
    geo->add_region(tmp);
  }

  return geo;
}

/// A 2-D pincell with slightly more complex division
static Geometry::SP_geometry
test_2D_pincell_via_factory(const size_t div, const size_t nrad)
{
  typedef RegionFactory RF;

  Point pitch(1.26);
  PinCell::vec_int mat_map(nrad+1, 0);
  PinCell::vec_dbl radii;
  if (nrad) radii = detran_utilities::linspace_center(0.0, 0.6, nrad);
  std::cout << radii.size() << std::endl;
  for (int i = 0; i < nrad; ++i) std::cout << radii[i] << std::endl;

  PinCell::SP_pincell pin;
  Point center(0, 0);
  pin = PinCell::Create(pitch, mat_map, radii, div, center);

  RF::vec_region regions;
  regions = RF::CreatePinCell(pin);

  Geometry::SP_geometry geo = Geometry::Create(1.26, 1.26, 0.0);
  for (size_t i = 0; i < regions.size(); ++i)
  {
    geo->add_region(regions[i]);
  }

  return geo;
}

/// A 2-D pincell assembly
static Geometry::SP_geometry test_2D_assembly(int n = 3)
{
  typedef RegionFactory RF;

  // Three cell types: moderator, fuel I, and fuel II
  Point P(1.26);
  PinCell::vec_int mat_M(1, 0);
  mat_M[0] = 0;
  PinCell::vec_int mat_I(3, 0);
  mat_I[0] = 1;
  mat_I[1] = 2;
  PinCell::vec_int mat_II(3, 0);
  mat_II[0] = 3;
  mat_II[1] = 4;
  PinCell::vec_dbl radii(2, 0);
  radii[0] = 0.50;
  radii[1] = 0.54;

  PinCell::SP_pincell M, F_I, F_II;
  M    = PinCell::Create(P, mat_M);
  F_I  = PinCell::Create(P, mat_I, radii, PinCell::DIVISION_HV_DIAG);
  F_II = PinCell::Create(P, mat_I, radii, PinCell::DIVISION_HV_DIAG);

  Assembly::SP_assembly a = Assembly::Create(n, n);
  int* pin_map_a;
  int pin_map_a_8[] =
  {
      2, 1, 1, 1, 1, 1, 1, 2,
      1, 0, 1, 1, 1, 1, 0, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 1, 1, 2, 2, 1, 1, 1,
      1, 1, 1, 2, 2, 1, 1, 1,
      1, 1, 1, 1, 1, 1, 1, 1,
      1, 0, 1, 1, 1, 1, 0, 1,
      2, 1, 1, 1, 1, 1, 1, 2
  };

  int pin_map_a_3[] =
  {
      1, 1, 1,
      1, 1, 1,
      1, 1, 1
  };
  if (n == 3)
    pin_map_a = pin_map_a_3;
  else
    pin_map_a = pin_map_a_8;

  Assembly::vec_int pin_map(n * n);
  for (int i = 0; i < pin_map.size(); ++i)
  {
    pin_map[i] = pin_map_a[i];
  }
  a->add_pincell(M);
  a->add_pincell(F_I);
  a->add_pincell(F_II);
  a->set_pincell_map(pin_map);

  RF::vec_region regions = RF::CreateAssembly(a);
  std::cout << " # reg = " << regions.size() << std::endl;

  Geometry::SP_geometry geo = Geometry::Create(n*1.26, n*1.26, 0.0);
  for (size_t i = 0; i < regions.size(); ++i)
  {
    geo->add_region(regions[i]);
  }

  return geo;
}

} // end namespace detran_geometry

#endif /* detran_geometry_CSG_FIXTURE_HH_ */

//----------------------------------------------------------------------------//
//              end of csg_fixture.cc
//----------------------------------------------------------------------------//

