/*
 * ReactionRates.cc
 *
 *  Created on: May 24, 2012
 *      Author: robertsj
 */

// Detran
#include "ReactionRates.hh"

// System
#include <algorithm>
#include <iostream>
#include <cmath>

namespace detran
{

// Constructor.
ReactionRates::ReactionRates(SP_material material,
                             SP_mesh mesh,
                             SP_state state)
  : b_material(material)
  , b_mesh(mesh)
  , b_state(state)
{
  /* ... */
}

vec_dbl ReactionRates::pin_power(double scale)
{

  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists("PINS"),
      "Mesh must have PINS map to compute pin power peaking factors.");

  return power("PINS", scale);

}

vec_dbl ReactionRates::assembly_power(double scale)
{

  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists("ASSEMBLIES"),
      "Mesh must have ASSEMBLIES map to compute pin power peaking factors.");

  return power("ASSEMBLIES", scale);

}

//---------------------------------------------------------------------------//
// Private Implementation
//---------------------------------------------------------------------------//

vec_dbl ReactionRates::power(std::string key, double scale)
{

  // We need a fission source.
  Insist(b_fissionsource,
      "Fission source is needed to compute powers.");

  // Compute the fission rate.
  vec_int mat_map = b_mesh->mesh_map("MATERIAL");
  State::moments_type fission_rate(b_mesh->number_cells(), 0.0);
  for(int g = 0; g < b_material->number_groups(); g++)
  {
    State::moments_type phi = b_state->phi(g);
    for (int cell = 0; cell < b_mesh->number_cells(); cell++)
    {
      fission_rate[cell] += phi[cell] *
                            b_material->sigma_f(mat_map[cell], g);
    }
  }

  // Get the region map of interest.
  vec_int map = b_mesh->mesh_map(key);

  // Number of unique regions.  This *assumes* sequential numbering,
  // which is true for the automated generation of pin cell and
  // assembly zones.
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
        int cell = b_mesh->index(i, j, k);
        double volume = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);
        power[map[cell]] += volume * fission_rate[cell];
        total_fission_rate += volume * fission_rate[cell];
      }
    }
  }

  // Adjust the scaling factor.
  scale /= total_fission_rate;

  // Scale the power and return.
  for (int region = 0; region < number; region++)
  {
    power[region] *= scale;
  }
  return power;
}

} // end namespace detran


