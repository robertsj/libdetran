//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGScatterFissionOperator.cc
 *  @brief MGScatterFissionOperator
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGScatterFissionOperator.hh"

namespace detran
{

//----------------------------------------------------------------------------//
MGScatterFissionOperator::
MGScatterFissionOperator(SP_input         input,
                         SP_material      material,
                         SP_mesh          mesh,
                         SP_scattersource S,
                         SP_fissionsource F,
                         size_t           cutoff,
                         bool             include_fission,
                         bool             adjoint)
  : Base(this)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_S(S)
  , d_F(F)
  , d_group_cutoff(cutoff)
  , d_include_fission(include_fission)
  , d_adjoint(adjoint)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  // Number of groups
  d_number_groups = d_material->number_groups();

  // Number of active groups
  int upper = d_number_groups;
  if (d_adjoint) upper = -1;
  d_number_active_groups = std::abs(upper - d_group_cutoff);

  d_moments_size = d_mesh->number_cells();


  // Set the operator size
  set_size(d_number_active_groups * d_moments_size);
}

//----------------------------------------------------------------------------//
void MGScatterFissionOperator::multiply(const Vector &V,  Vector &V_out)
{
  // \todo This doesn't look like it is adjoint-ready

  size_t size_moments = d_mesh->number_cells();

  // Copy input vector to a multigroup flux; only the Krylov block is used.
  State::vec_moments_type phi(d_number_groups,
      State::moments_type(size_moments, 0.0));
  for (int g = d_group_cutoff; g < d_number_groups; ++g)
    for (int i = 0; i < size_moments; ++i)
      phi[g][i] = V[(g - d_group_cutoff) * size_moments + i];

  // Construct the action of the scattering and fission operators
  for (int g = d_group_cutoff; g < d_number_groups; ++g)
  {
    State::moments_type source(size_moments, 0.0);
    // Add scatter
    d_S->build_total_group_source(g, d_group_cutoff, phi, source);
    // Add fission
    if (d_include_fission)
      d_F->build_total_group_source(g, phi, source);
    for (int i = 0; i < size_moments; i++)
      V_out[(g - d_group_cutoff) * size_moments + i] = source[i];
  }

}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGScatterFissionOperator.cc
//----------------------------------------------------------------------------//




