//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ReactionRates.cc
 *  @author robertsj
 *  @date   May 24, 2012
 *  @brief  ReactionRates member definitions
 */
//---------------------------------------------------------------------------//

#include "ReactionRates.hh"
#include <algorithm>
#include <iostream>
#include <cmath>

namespace detran_postprocess
{

//---------------------------------------------------------------------------//
ReactionRates::ReactionRates(SP_material material,
                             SP_mesh mesh,
                             SP_state state)
  : b_material(material)
  , b_mesh(mesh)
  , b_state(state)
{
  /* ... */
}

//---------------------------------------------------------------------------//
detran_utilities::vec_dbl
ReactionRates::region_power(std::string key, double scale)
{

  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists(key),
      "The mesh map " + key + " does not exist.");

  // Compute the fission rate.
  vec_int mat_map = b_mesh->mesh_map("MATERIAL");
  detran::State::moments_type fission_rate(b_mesh->number_cells(), 0.0);
  for(size_t g = 0; g < b_material->number_groups(); g++)
  {
    detran::State::moments_type phi = b_state->phi(g);
    for (int cell = 0; cell < b_mesh->number_cells(); cell++)
    {
      fission_rate[cell] += phi[cell] *
                            b_material->sigma_f(mat_map[cell], g);
    }
  }

  // Get the region map of interest.
  vec_int map = b_mesh->mesh_map(key);

  // Get the maximum index, which we take to the be one less
  // then the number of unique regions.  If there is non-sequential
  // number, then elements of the power vector are simply zero.
  // Number of unique regions.  This *assumes* sequential numbering,
  int number = 1 + *std::max_element(map.begin(), map.end());

  // Compute the (unscaled) region powers.
  vec_dbl power(number, 0.0);
  double total_fission_rate = 0.0;
  for (int k = 0; k < b_mesh->number_cells_z(); k++)
  {
    for (int j = 0; j < b_mesh->number_cells_y(); j++)
    {
      for (int i = 0; i < b_mesh->number_cells_x(); i++)
      {
        // Define the mesh cell and volume
        int cell = b_mesh->index(i, j, k);
        double volume = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);

        // Add contribution to this region power
        power[map[cell]] += volume * fission_rate[cell];

        // Add contribution to the total power
        total_fission_rate += volume * fission_rate[cell];
      }
    }
  }

  if (scale > 0)
  {
    // Adjust the scaling factor.
    scale /= total_fission_rate;

    // Scale the power and return.
    for (int region = 0; region < number; region++)
    {
      power[region] *= scale;
    }
  }

  return power;
}

//---------------------------------------------------------------------------//
detran_utilities::vec_dbl
ReactionRates::edit(std::string key,
                    const vec_dbl &fine_mesh_function,
                    bool mean)
{
  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists(key),
      "The mesh map " + key + " does not exist.");

  // Get the region map of interest.
  vec_int map = b_mesh->mesh_map(key);

  // Number of unique regions.  This *assumes* sequential numbering,
  int number = 1 + *std::max_element(map.begin(), map.end());

  vec_dbl edit_region_function(number, 0.0);
  vec_dbl edit_region_volume(number, 0.0);

  for (size_t k = 0; k < b_mesh->number_cells_z(); k++)
  {
    for (size_t j = 0; j < b_mesh->number_cells_y(); j++)
    {
      for (size_t i = 0; i < b_mesh->number_cells_x(); i++)
      {
        // Define the mesh cell and volume
        int cell = b_mesh->index(i, j, k);
        double volume = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);

        // Add contribution to this edit region
        edit_region_function[map[cell]] += volume * fine_mesh_function[cell];
        edit_region_volume[map[cell]] += volume;
      }
    }
  }

  if (mean)
    for (int cell = 0; cell < number; ++cell)
      edit_region_function[map[cell]] /= edit_region_volume[map[cell]];

  return edit_region_function;
}


} // end namespace detran_postprocess

//---------------------------------------------------------------------------//
//              end of ReactionRates.cc
//---------------------------------------------------------------------------//
