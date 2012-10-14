//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   SweepSource.i.hh
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  SweepSource inline member definitions.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPSOURCE_I_HH_
#define SWEEPSOURCE_I_HH_

#include <iostream>
#include <omp.h>

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::build_fixed(const size_t g)
{
  // Zero out moments source.
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    d_fixed_group_source[cell] = 0.0;
  }
  // Add external sources, if present.
  for (int i = 0; i < d_moment_external_sources.size(); i++)
  {
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      d_fixed_group_source[cell] +=
        d_moment_external_sources[i]->source(cell, g);
    }
  }
  // Add fission source if present
  if (d_fissionsource and !d_implicit_fission)
  {
    const State::moments_type &qf = d_fissionsource->source(g);
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      d_fixed_group_source[cell] += qf[cell];
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
build_fixed_with_scatter(const size_t g)
{
  // Add the external and/or fission source first.
  build_fixed(g);
  // Add the in-scatter.
  d_scattersource->build_in_scatter_source(g, d_fixed_group_source);
  // Add the in-fission if applicable.
  if (d_implicit_fission)
    d_fissionsource->build_in_fission_source(g, d_fixed_group_source);
}

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
build_fixed_with_downscatter(const size_t g, const size_t g_cutoff)
{
  // Add the external and/or fission source first.
  build_fixed(g);
  // Add the in-scatter.
  d_scattersource->build_downscatter_source(g, g_cutoff, d_fixed_group_source);
  // Note, we're not adding fission here, since a multiplying problem
  // should never be solved in downscatter mode
}

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
build_within_group_scatter(const size_t g, const moments_type &phi)
{
  // Zero out moments source.
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    d_scatter_group_source[cell] = 0.0;
  // Build within-group scattering
  d_scattersource->build_within_group_source(g, phi, d_scatter_group_source);
  if (d_implicit_fission)
    d_fissionsource->build_within_group_source(g, phi, d_scatter_group_source);
}

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
build_total_scatter(const size_t g, const size_t g_cutoff,
                    const State::vec_moments_type &phi)
{
  // Zero out moments source.
  d_scatter_group_source.assign(phi[g].size(), 0.0);
  // Build total scattering source
  d_scattersource->build_total_group_source(g, g_cutoff, phi,
                                            d_scatter_group_source);
  if (d_implicit_fission)
  {
    // For multiplying problems, there should be no downscatter block.
    Assert(g_cutoff == 0);
    d_fissionsource->build_total_group_source(g, phi,
                                              d_scatter_group_source);
  }
}

//---------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
source(const size_t g, const size_t o, const size_t a, sweep_source_type &s)
{
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    s[cell] = 0.0;

  // Add moment contributions.
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    // Add fixed contribution
    s[cell] += d_fixed_group_source[cell] *
                      (*d_MtoD)(o, a, 0, 0);
    // Add scatter contribution
    s[cell] += d_scatter_group_source[cell] *
                      (*d_MtoD)(o, a, 0, 0);
  }

  // Add discrete contributions if present.
  for (int i = 0; i < d_discrete_external_sources.size(); i++)
  {
    for (int cell = 0; cell < d_mesh->number_cells(); cell++)
    {
      s[cell] += d_discrete_external_sources[i]->source(cell, g);
    }
  }

}

//---------------------------------------------------------------------------//
template <class D>
void SweepSource<D>::reset()
{
  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
  {
    d_fixed_group_source[cell] = 0.0;
    d_scatter_group_source[cell] = 0.0;
  }
}

// Explicit instantiations
template class SweepSource<_1D>;
template class SweepSource<_2D>;
template class SweepSource<_3D>;

} // namespace detran

#endif /* SWEEPSOURCE_I_HH_ */

//---------------------------------------------------------------------------//
//              end of SweepSource.i.hh
//---------------------------------------------------------------------------//
