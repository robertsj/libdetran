//----------------------------------*-C++-*----------------------------------//
/**
 * @file   FissionSource.cc
 * @brief  FissionSource class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
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
  , d_adjoint(false)
{
  Require(d_state);
  Require(d_mesh);
  Require(d_material);

  d_mat_map = d_mesh->mesh_map("MATERIAL");
  d_number_groups = d_material->number_groups();
  d_density.assign(d_mesh->number_cells(), 0.0);
  d_source.resize(d_material->number_groups(),
                  moments_type(mesh->number_cells(), 0.0));
  if (d_state->get_input()->check("adjoint"))
    d_adjoint = 0 != d_state->get_input()->get<int>("adjoint");
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
  int ng = d_material->number_groups();
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    d_density[cell] = 0.0;
    for (int g = 0; g < ng; g++)
    {
      d_density[cell] += d_material->nu_sigma_f(d_mat_map[cell], g);
    }
  }
  double norm_density = detran_utilities::norm(d_density, "L1");
  Ensure(norm_density > 0.0);
  detran_utilities::vec_scale(d_density, 1.0/norm_density);
}

//---------------------------------------------------------------------------//
void FissionSource::set_scale(const double scale)
{
  d_scale = scale;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of FissionSource.cc
//---------------------------------------------------------------------------//
