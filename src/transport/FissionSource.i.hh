//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  FissionSource.i.hh
 *  @brief FissionSource inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_FISSIONSOURCE_I_HH_
#define detran_FISSIONSOURCE_I_HH_

namespace detran
{

//----------------------------------------------------------------------------//
inline void FissionSource::setup_outer(const double scale)
{
  d_scale = scale;
  for (size_t g = 0; g < d_material->number_groups(); ++g)
  {
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      double v = d_adjoint ? d_material->nu_sigma_f(d_mat_map[cell], g)
                           : d_material->chi(d_mat_map[cell], g);
      d_source[g][cell] = d_scale * d_density[cell] * v;
    }
  }
}


//----------------------------------------------------------------------------//
inline void FissionSource::update()
{
  d_density.assign(d_density.size(), 0.0);
  for(size_t g = 0; g < d_number_groups; g++)
  {
    State::moments_type phi = d_state->phi(g);
    for (size_t cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      double v = d_adjoint ? d_material->chi(d_mat_map[cell], g)
                           : d_material->nu_sigma_f(d_mat_map[cell], g);
      d_density[cell] += phi[cell] * v;
    }
  }
}

//----------------------------------------------------------------------------//
inline const State::moments_type& FissionSource::source(const size_t g) const
{
  Require(g < d_number_groups);
  return d_source[g];
}

//----------------------------------------------------------------------------//
inline const State::moments_type& FissionSource::density() const
{
  return d_density;
}

//----------------------------------------------------------------------------//
inline State::moments_type& FissionSource::density()
{
  return d_density;
}


//----------------------------------------------------------------------------//
// INTERFACE MIMICING SCATTER SOURCE
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
inline void FissionSource::
build_within_group_source(const size_t               g,
                          const State::moments_type &phi,
                          State::moments_type       &source)
{
  Require(g < d_material->number_groups());
  Require(phi.size() == source.size());

  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    source[cell] += phi[cell] * d_scale *
                    d_material->chi(d_mat_map[cell], g) *
                    d_material->nu_sigma_f(d_mat_map[cell], g);
  }
}

//----------------------------------------------------------------------------//
inline void FissionSource::
build_in_fission_source(const size_t  g,
                        moments_type &source)
{
  Require(g < d_material->number_groups());

  for (size_t gp = 0; gp < d_material->number_groups(); ++gp)
  {
    if (gp == g) continue;
    const size_t g_f = g_from(g, gp);
    const size_t g_t = g_to(g, gp);
    const moments_type &phi = d_state->phi(gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      source[cell] += phi[cell] * d_scale *
                      d_material->chi(d_mat_map[cell], g_t) *
                      d_material->nu_sigma_f(d_mat_map[cell], g_f);
    }
  }
}

//----------------------------------------------------------------------------//
inline void FissionSource::
build_total_group_source(const size_t                   g,
                         const State::vec_moments_type &phi,
                         State::moments_type           &source)
{
  Require(g < d_material->number_groups());

  for (size_t gp = 0; gp < d_material->number_groups(); ++gp)
  {
    const size_t g_f = g_from(g, gp);
    const size_t g_t = g_to(g, gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      source[cell] += phi[gp][cell] * d_scale *
                      d_material->chi(d_mat_map[cell], g_t) *
                      d_material->nu_sigma_f(d_mat_map[cell], g_f);
    }
  }
}

//----------------------------------------------------------------------------//
inline FissionSource::size_t
FissionSource::g_from(const size_t g, const size_t gp) const
{
  return d_adjoint ? g : gp;
}

//----------------------------------------------------------------------------//
inline FissionSource::size_t
FissionSource::g_to(const size_t g, const size_t gp) const
{
  return d_adjoint ? gp : g;
}

} // namespace detran

#endif /* detran_FISSIONSOURCE_I_HH_ */

//----------------------------------------------------------------------------//
//              end of FissionSource.i.hh
//----------------------------------------------------------------------------//
