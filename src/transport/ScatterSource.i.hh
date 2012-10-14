//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ScatterSource.i.hh
 *  @author robertsj
 *  @date   Sep 7, 2012
 *  @brief  ScatterSource inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef detran_SCATTERSOURCE_I_HH_
#define detran_SCATTERSOURCE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline void ScatterSource::
build_within_group_source(const size_t g,
                          const moments_type &phi,
                          moments_type &source)
{
  // Preconditions
  Require(g < d_material->number_groups());
  Require(phi.size() == source.size());

  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, g);
  }
}

//---------------------------------------------------------------------------//
inline void ScatterSource::
build_in_scatter_source(const size_t g,
                        moments_type &source)
{
  // Preconditions.
  Require(g < d_material->number_groups());

  // Add downscatter.
  for (int gp = d_material->lower(g); gp < g; gp++) //
  {
    const moments_type &phi = d_state->phi(gp);
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
    }
  }
  // Add upscatter.
  for (int gp = g + 1; gp <= d_material->upper(g); gp++)
  {
    const moments_type &phi = d_state->phi(gp);
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
    }
  }
}

//---------------------------------------------------------------------------//
inline void ScatterSource::
build_downscatter_source(const size_t g,
                         const size_t g_cutoff,
                         moments_type &source)
{
  // Preconditions.
  Require(g < d_material->number_groups());
  Require(g_cutoff <= d_material->number_groups());

  // Add downscatter.
  for (int gp = d_material->lower(g); gp < g_cutoff; gp++) //
  {
    moments_type phi = d_state->phi(gp);
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      source[cell] += phi[cell] * d_material->sigma_s(d_mat_map[cell], g, gp);
    }
  }

}

//---------------------------------------------------------------------------//
inline void ScatterSource::
build_total_group_source(const size_t g,
                         const size_t g_cutoff,
                         const State::vec_moments_type &phi,
                         moments_type &source)
{
  // Preconditions
  Require(g < d_material->number_groups());

  for (int gp = g_cutoff; gp <= d_material->upper(g); gp++)
  {
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      source[cell] += phi[gp][cell] *
          d_material->sigma_s(d_mat_map[cell], g, gp);
    }
  }
}

} // end namespace detran

#endif /* detran_SCATTERSOURCE_I_HH_ */
