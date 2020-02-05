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
                         size_t           sf_switch ,
                         bool             adjoint)
  : Base(this)
  , d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_S(S)
  , d_F(F)
  , d_group_cutoff(cutoff)
  , d_include_scatter(sf_switch < 2)
  , d_include_fission(sf_switch > 0)
  , d_adjoint(adjoint)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);

  using detran_utilities::range;

  // Number of groups
  d_number_groups = d_material->number_groups();

  // Number of active groups
  size_t lower = d_group_cutoff;
  size_t upper = d_number_groups;
  if (d_adjoint) upper = -1;

  d_number_active_groups = std::abs((int)upper - (int)d_group_cutoff);

  d_groups = range<size_t>(lower, upper);

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

  groups_t::iterator g = d_groups.begin();
  for (; g != d_groups.end(); ++g)
    for (int i = 0; i < size_moments; ++i)
      phi[*g][i] = V[(*g - d_group_cutoff) * size_moments + i];

  // Construct the action of the scattering and/or fission operators
  g = d_groups.begin();
  for (; g != d_groups.end(); ++g)
  {
    State::moments_type source(size_moments, 0.0);
    // Add scatter
    if (d_include_scatter)
      d_S->build_total_group_source(*g, d_group_cutoff, phi, source);
    // Add fission
    if (d_include_fission)
      d_F->build_total_group_source(*g, phi, source);
    for (int i = 0; i < size_moments; ++i)
      V_out[(*g - d_group_cutoff) * size_moments + i] = source[i];
  }
}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGScatterFissionOperator.cc
//----------------------------------------------------------------------------//
