//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PulsedExternalSource.cc
 *  @brief  PulsedExternalSource
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#include "PulsedExternalSource.hh"

namespace detran
{

//---------------------------------------------------------------------------//
PulsedExternalSource::PulsedExternalSource(const size_t       number_groups,
                                           SP_mesh            mesh,
                                           SP_externalsource  fixed_source,
                                           const double       peak_time,
                                           const double       fwhm,
                                           bool               discrete)
  : Base(number_groups, mesh, SP_quadrature(0), discrete)
  , d_fixed_source(fixed_source)
  , d_peak_time(peak_time)
  , d_fwhm(fwhm)
  , d_factor(0.0)
{
  // Preconditions
  Require(d_fixed_source);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PulsedExternalSource.cc
//---------------------------------------------------------------------------//
