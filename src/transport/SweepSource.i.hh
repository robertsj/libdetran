//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SweepSource.i.hh
 *  @brief SweepSource inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_SWEEPSOURCE_I_HH_
#define detran_SWEEPSOURCE_I_HH_

#include <iostream>
#ifdef DETRAN_ENABLE_OPENMP
// Xcode apparently doesn't have OpenMP by default
#include <omp.h>
#endif

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::build_fixed(const size_t g)
{
  // Zero out moments source.
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    d_fixed_group_source[cell] = 0.0;
  }
  // Add external sources, if present.
  for (size_t i = 0; i < d_moment_external_sources.size(); ++i)
  {
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      d_fixed_group_source[cell] +=
        d_moment_external_sources[i]->source(cell, g);
    }
  }
  // Add fission source if present
  if (d_fissionsource && !d_implicit_fission)
  {
    const State::moments_type &qf = d_fissionsource->source(g);
    for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    {
      d_fixed_group_source[cell] += qf[cell];
    }
  }
}

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
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

//----------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
build_within_group_scatter(const size_t g, const moments_type &phi)
{
  // Zero out moments source.
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    d_scatter_group_source[cell] = 0.0;
  // Build within-group scattering
  d_scattersource->build_within_group_source(g, phi, d_scatter_group_source);
  if (d_implicit_fission)
    d_fissionsource->build_within_group_source(g, phi, d_scatter_group_source);
}

//----------------------------------------------------------------------------//
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
    Assert(g_cutoff == 0 || g_cutoff == d_state->number_groups() - 1);
    d_fissionsource->build_total_group_source(g, phi,
                                              d_scatter_group_source);
  }
}

//----------------------------------------------------------------------------//
template <class D>
inline void SweepSource<D>::
source(const size_t g, const size_t o, const size_t a, sweep_source_type &s)
{
  // \todo A study on the optimal implementation is warranted.  This
  //       will certainly be group/quad dependent.

  double *s_a = &s[0];
  double *fixed_a = &d_fixed_group_source[0];
  double *scatter_a = &d_scatter_group_source[0];
  double mtod = (*d_MtoD)(o, a, 0, 0);

//  // Add fixed contributions
//  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
//    s_a[cell] = fixed_a[cell] * mtod;
//
//  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
//    s_a[cell] += scatter_a[cell] * mtod;

  // Add fixed contributions
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
    s_a[cell] = (fixed_a[cell] + scatter_a[cell])* mtod;
//
//  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
//    s[cell] = 0.0;
//
//
//  // Add moment contributions.
//  for (int cell = 0; cell < d_mesh->number_cells(); cell++)
//  {
//    // Add fixed contribution
//    s[cell] += d_fixed_group_source[cell] *
//                      (*d_MtoD)(o, a, 0, 0);
//    // Add scatter contribution
//    s[cell] += d_scatter_group_source[cell] *
//                      (*d_MtoD)(o, a, 0, 0);
//  }

  // Add discrete contributions if present.
  if (d_discrete_external_source_flag)
  {
    size_t angle = d_quadrature->index(o, a);
    for (size_t i = 0; i < d_discrete_external_sources.size(); ++i)
      for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
        s[cell] += d_discrete_external_sources[i]->source(cell, g, angle);
  }

}

//----------------------------------------------------------------------------//
template <class D>
void SweepSource<D>::reset()
{
  for (size_t cell = 0; cell < d_mesh->number_cells(); ++cell)
  {
    d_fixed_group_source[cell] = 0.0;
    d_scatter_group_source[cell] = 0.0;
  }
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

TRANSPORT_INSTANTIATE_EXPORT(SweepSource<_1D>)
TRANSPORT_INSTANTIATE_EXPORT(SweepSource<_2D>)
TRANSPORT_INSTANTIATE_EXPORT(SweepSource<_3D>)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<SweepSource<_1D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<SweepSource<_2D> >)
TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<SweepSource<_3D> >)

} // namespace detran

#endif /* detran_SWEEPSOURCE_I_HH_ */

//----------------------------------------------------------------------------//
//              end of SweepSource.i.hh
//----------------------------------------------------------------------------//
