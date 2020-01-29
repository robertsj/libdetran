//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PinCell.hh
 *  @brief PinCell class member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts  
 */
//----------------------------------------------------------------------------//

#include "PinCell.hh"
#include "utilities/SoftEquivalence.hh"
#include "utilities/Warning.hh"
#include <cmath>
#include <iostream>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
PinCell::PinCell(const Point      &pitch,
                 const vec_int    &mat_map,
                 const vec_dbl    &radii,
                 const size_t      division,
                 const Point      &pincenter)
  : d_pitch(pitch)
  , d_mat_map(mat_map)
  , d_division(division)
  , d_radii(radii)
  , d_pin_center(pincenter)
{
  Insist(d_pitch.x() > 0.0, "Pitch must be positive");
  Insist(d_mat_map.size() == d_radii.size() + 1,
         "Material id vector size != to # of pin + moderator regions");
  Insist(d_division < END_DIVISION_SCHEMES, "Invalid division scheme");
  if (d_radii.size())
  {
    Insist(d_radii[0] > 0.0, "Radii must be positive");
  }
  std::cout << " # radii = " << radii.size() << std::endl;
  for (size_t i = 1; i < d_radii.size(); ++i)
  {
    Insist(d_radii[i] > d_radii[i-1], "Radii must be monotonic increasing:"
           + AsString(d_radii[i-1]) + " followed by " + AsString(d_radii[i]));
  }
  if (d_pitch.y() <= 0.0)
  {
    d_pitch = Point(d_pitch.x(), d_pitch.x());
  }
}

//----------------------------------------------------------------------------//
PinCell::SP_pincell PinCell::Create(const Point      &pitch,
                                    const vec_int    &mat_map,
                                    const vec_dbl    &radii,
                                    const size_t      division,
                                    const Point      &pincenter)
{
  SP_pincell p(new PinCell(pitch, mat_map, radii, division, pincenter));
  return p;
}

//----------------------------------------------------------------------------//
void PinCell::meshify(int number_meshes, bool flag)
{
  Require(number_meshes > 0);
  if (flag && d_radii.size() != 1)
  {
    // Best fit mesh only for single radius.
    detran_utilities::warning(detran_utilities::USER_INPUT,
      "Best fit pincell only for a single radius.");
    flag = false;
  }
  if (flag && number_meshes < 7)
  {
    // Best fit mesh only for single radius.
    detran_utilities::warning(detran_utilities::USER_INPUT,
      "Best fit pincell needs at least 7 meshes.");
    flag = false;
  }

  // Use cell-center material and uniform mesh
  if (!flag)
  {

    // Fine mesh width.
    double width = d_pitch.x() / number_meshes;

    // Fine mesh edges.  Assumes constant spacing.
    vec_dbl edges(number_meshes + 1, 0.0);

    // Number of cells in the meshed pin cell.
    int number_cells = number_meshes*number_meshes;

    // Temporary fine mesh material map
    vec_int tmp_mat_map(number_cells, 1);
    vec_int tmp_reg_map(number_cells, 1);

    //double volume = d_number_regions, 1);
    for (int j = 0; j < number_meshes; j++)
    {
      edges[j+1] = edges[j] + width;
      for (int i = 0; i < number_meshes; i++)
      {
        // Which region am I in?
        int r = find_region(i, j, width);
        // Which material does this region have?
        int m = d_mat_map[r];
        // Assign the values.
        int cell = i + j * number_meshes;
        tmp_mat_map[cell] = m;
        tmp_reg_map[cell] = r;
      }
    }

    // Create my mesh.
    d_mesh = new Mesh2D(edges, edges, tmp_mat_map);

    // Add maps
    d_mesh->add_mesh_map("MATERIAL", tmp_mat_map);
    d_mesh->add_mesh_map("REGION", tmp_reg_map);

  }
  // Use volume-conserving approximate cylinder
  else
  {
    // Best fit parameters. 1.604968750000000   0.176223981031790
    double L = 1.60496875 * d_radii[0];
    double delta = 0.176223981031790 * d_radii[0];
    // Dimensions
    double half_pitch = 0.5*d_pitch.x();
    double width_mod = half_pitch-0.5*L-delta;

    // Number of fine meshes in moderator, delta, and square
    // regions, along the x axis
    int num_mod    =
      std::max(std::floor(0.5*number_meshes * width_mod / half_pitch), 1.0);
    int num_delta  =
      std::max(std::floor(0.5 * number_meshes * delta / half_pitch),   1.0);
    int num_L_4    =
      std::max(std::floor(0.5 * number_meshes * 0.25*L / half_pitch), 1.0);
    int num_L_2    =
      std::max(std::floor(0.25 * number_meshes * L / half_pitch)-1.0,
               double(number_meshes-2.0*(num_mod+num_delta+num_L_4)));

//    std::cout << " num_mod = "      << num_mod      << std::endl;
//    std::cout << " num_delta = "    << num_delta    << std::endl;
//    std::cout << " num_L_4 = "      << num_L_4      << std::endl;
//    std::cout << " num_L_2 = "      << num_L_2      << std::endl;

    int num_mesh = 2*num_mod + 2*num_delta + 2*num_L_4 + num_L_2;

    if ( (number_meshes-num_mesh) > 0)
    {
      num_L_2 += 1;
      num_mesh += 1;
    }
    if ( (number_meshes-num_mesh) > 1)
    {
      num_L_4 += 1;
      num_mesh += 2;
    }
    if ( (number_meshes-num_mesh) > 1)
    {
      num_mod += 1;
      num_mesh += 2;
    }
    if ( (number_meshes-num_mesh) > 1)
    {
      num_delta += 1;
      num_mesh += 2;
    }
    if ( (number_meshes-num_mesh) > 0)
    {
      num_L_2 += 1;
      num_mesh += 1;
    }

    Ensure(num_mesh == number_meshes);
    Ensure(num_mod > 0);
    Ensure(num_delta > 0);
    Ensure(num_L_4 > 0);
    Ensure(num_L_2 > 0);

    // Coarse mesh edges
    vec_dbl coarse(8, 0.0);
    // Fine mesh counts
    vec_int fine(7, 0);
    // Coarse mesh materials and region
    vec_int material(49, d_mat_map[1]);
    vec_int regions(49, 1);

    // Set the fine mesh counts for the 5x5 coarse map.
    fine[0] = num_mod;
    fine[1] = num_delta;
    fine[2] = num_L_4;
    fine[3] = num_L_2;
    fine[4] = num_L_4;
    fine[5] = num_delta;
    fine[6] = num_mod;

    // Set the edges.
    coarse[0] = 0.0;
    coarse[1] = width_mod;
    coarse[2] = coarse[1] + delta;
    coarse[3] = coarse[2] + 0.25*L;
    coarse[4] = coarse[3] + 0.5*L;
    coarse[5] = coarse[4] + 0.25*L;
    coarse[6] = coarse[5] + delta;
    coarse[7] = coarse[6] + width_mod;
    Ensure(detran_utilities::soft_equiv(coarse[7], d_pitch.x()));

    // Set the material.  Default is zero, which
    // is fuel.  We need to set the 13 fuel cells.
    int fuel_index[] = {10,16,17,18,22,23,24,25,26,30,31,32,38};
    for (int i = 0; i < 13; i++)
    {
      material[fuel_index[i]] = d_mat_map[0];
      regions[fuel_index[i]] = 0;
    }

    // Create my mesh.
    d_mesh = new Mesh2D(fine, fine, coarse, coarse, material);
    d_mesh->add_coarse_mesh_map("REGION", regions);

  }

}

//----------------------------------------------------------------------------//
int PinCell::find_region(int i, int j, double width)
{
  double x = (i + 0.5) * width;
  double y = (j + 0.5) * width;
  // Loop through the radii. If I'm in there, that's where I live.
  int r = d_radii.size(); // start off in outer region.
  double hp = 0.5 * d_pitch.x();
  for (size_t p = 0; p < d_radii.size(); p++)
  {
    if (std::sqrt((x - hp) * (x - hp) + (y - hp) * (y - hp)) < d_radii[p])
    {
      r = p;
      break;
    }
  }
  return r;
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of PinCell.cc
//----------------------------------------------------------------------------//

