//----------------------------------*-C++-*----------------------------------//
/**
 * @file   FissionSource.cc
 * @author robertsj
 * @date   Apr 10, 2012
 * @brief  FissionSource class definition.
 */
//---------------------------------------------------------------------------//

#include "transport/FissionSource.hh"
#include "utilities/MathUtilities.hh"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
FissionSource::FissionSource(SP_state state,
                             SP_mesh mesh,
                             SP_material material)
  : d_state(state)
  , d_mesh(mesh)
  , d_material(material)
  , d_scale(1.0)
{
  // Preconditions
  Require(d_state);
  Require(d_mesh);
  Require(d_material);

  d_number_groups = d_material->number_groups();
  d_density.assign(d_mesh->number_cells(), 0.0);
  d_source.resize(d_material->number_groups(),
                  moments_type(mesh->number_cells(), 0.0));
}

//---------------------------------------------------------------------------//
FissionSource::SP_fissionsource
FissionSource::Create(SP_state state, SP_mesh mesh, SP_material material)
{
  SP_fissionsource p;
  p = new FissionSource(state, mesh, material);
  return p;
}

//---------------------------------------------------------------------------//
void FissionSource::initialize()
{
  /*
   *  Define the density as the normalized sum of
   *  nu * fission cross section.  Normalized
   *  using the L1 norm.
   */
  vec_int mat_map = d_mesh->mesh_map("MATERIAL");
  int ng = d_material->number_groups();
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    d_density[cell] = 0.0;
    for (int g = 0; g < ng; g++)
    {
      d_density[cell] += d_material->nu_sigma_f(mat_map[cell], g);
    }
  }
  double norm_density = detran_utilities::norm(d_density, "L1");
  Require(norm_density > 0.0);
  detran_utilities::vec_scale(d_density, 1.0/norm_density);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of FissionSource.cc
//---------------------------------------------------------------------------//
