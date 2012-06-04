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

vec_dbl ReactionRates::pin_peaking()
{

  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists("PINS"),
      "Mesh must have PINS map to compute pin power peaking factors.");

  return peaking("PINS");

}

vec_dbl ReactionRates::assembly_peaking()
{

  // Make sure the mesh has what we need.
  Insist(b_mesh->mesh_map_exists("ASSEMBLIES"),
      "Mesh must have ASSEMBLIES map to compute pin power peaking factors.");

  return peaking("ASSEMBLIES");

}


//---------------------------------------------------------------------------//
// Private Implementation
//---------------------------------------------------------------------------//

vec_dbl ReactionRates::peaking(std::string key)
{

  // We need a fission source.
  Insist(b_fissionsource,
      "Fission source is needed to compute pin power peaking factors.");

  // Get the pin or assembly map.
  vec_int map = b_mesh->mesh_map(key);

  std::cout << " number = "  << *std::max_element(map.begin(), map.end()) << std::endl;
  int number = 1 + *std::max_element(map.begin(), map.end());

  // Compute pin power peaking factors
  State::moments_type fission_rate = b_fissionsource->density();
  double total_fission_rate = 0;
  // std::accumulate(fission_rate.begin(), fission_rate.end(), 0.0);

  vec_dbl pf(number, 0.0);

  for (int k = 0; k < b_mesh->number_cells_z(); k++)
  {
    for (int j = 0; j < b_mesh->number_cells_y(); j++)
    {
      for (int i = 0; i < b_mesh->number_cells_x(); i++)
      {
        int cell = b_mesh->index(i, j, k);
        double V = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);
        pf[map[cell]] += V * fission_rate[cell];
        total_fission_rate += V * fission_rate[cell];
      }
    }
  }

  int fissile = 0;
  for (int p = 0; p < number; p++)
  {
    // Don't count pins and assemblies without a (real) fission rate
    if (pf[p] > 0.00001) fissile++;
  }
  double average = total_fission_rate / double(fissile);

  for (int p = 0; p < number; p++)
  {
    pf[p] /= average;
  }

  // Returning pin power so total power = number of fissile pins
  return pf;
}



} // end namespace detran


