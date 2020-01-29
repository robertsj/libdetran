//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ScatterSource.i.hh
 *  @brief ScatterSource inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SCATTERSOURCE_I_HH_
#define detran_SCATTERSOURCE_I_HH_

namespace detran
{

//----------------------------------------------------------------------------//
inline void ScatterSource::
build_within_group_source(const size_t        g,
                          const moments_type &phi,
                          moments_type       &s)
{
  Require(g < d_material->number_groups());
  Require(phi.size() == s.size());

  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    s[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, g);
  }
}

//----------------------------------------------------------------------------//
inline void ScatterSource::
build_in_scatter_source(const size_t  g,
                        moments_type &s)
{
  Require(g < d_material->number_groups());

  groups_iter gp = groups(lower(g), upper(g), true);
  for (; gp != d_groups.end(); ++gp)
  {
    if (g == *gp) continue;
    const size_t g_f = g_from(g, *gp);
    const size_t g_t = g_to(g, *gp);
    const moments_type &phi = d_state->phi(*gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      s[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g_t, g_f);
    }
  }
}

//----------------------------------------------------------------------------//
inline void ScatterSource::
build_downscatter_source(const size_t  g,
                         const size_t  g_cutoff,
                         moments_type &s)
{
  Require(g < d_material->number_groups());
  Require(g_cutoff <= d_material->number_groups());

  groups_iter gp = groups(lower(g), g_cutoff, false);
  for (; gp != d_groups.end(); ++gp)
  {
    const size_t g_f = g_from(g, *gp);
    const size_t g_t = g_to(g, *gp);
    moments_type &phi = d_state->phi(*gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      s[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g_t, g_f);
    }
  }
}

//----------------------------------------------------------------------------//
inline void ScatterSource::
build_total_group_source(const size_t                   g,
                         const size_t                   g_cutoff,
                         const State::vec_moments_type &phi,
                         moments_type                  &s)
{
  Require(g < d_material->number_groups());

  groups_iter gp = groups(g_cutoff, upper(g), true);
  for (; gp != d_groups.end(); ++gp)
  {
    const size_t g_f = g_from(g, *gp);
    const size_t g_t = g_to(g, *gp);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      s[cell] += phi[*gp][cell] *
                 d_material->sigma_s(d_mat_map[cell], g_t, g_f);
    }
  }
}

//----------------------------------------------------------------------------//
inline ScatterSource::size_t ScatterSource::lower(const size_t g) const
{
  return d_material->lower(g, d_adjoint);
}

//----------------------------------------------------------------------------//
inline ScatterSource::size_t ScatterSource::upper(const size_t g) const
{
  return d_material->upper(g, d_adjoint);
}

//----------------------------------------------------------------------------//
inline ScatterSource::size_t
ScatterSource::g_from(const size_t g, const size_t gp) const
{
  return d_adjoint ? g : gp;
}

//----------------------------------------------------------------------------//
inline ScatterSource::size_t
ScatterSource::g_to(const size_t g, const size_t gp) const
{
  return d_adjoint ? gp : g;
}

//----------------------------------------------------------------------------//
inline ScatterSource::groups_iter
ScatterSource::groups(const size_t g_start, const size_t g_finish, bool inc)
{
  d_groups = detran_utilities::range<size_t>(g_start, g_finish, inc);
  return d_groups.begin();
}

} // end namespace detran

#endif /* detran_SCATTERSOURCE_I_HH_ */

//----------------------------------------------------------------------------//
//              end of ScatterSource.i.hh
//----------------------------------------------------------------------------//
