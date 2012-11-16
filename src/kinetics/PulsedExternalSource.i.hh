//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PulsedExternalSource.i.hh
 *  @brief  PulsedExternalSource.i
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_PULSEDEXTERNALSOURCE_I_HH_
#define detran_PULSEDEXTERNALSOURCE_I_HH_

#include <cmath>
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
inline double PulsedExternalSource::source(const size_t cell,
                                           const size_t group)
{
  // Preconditions
  Require(cell  < d_mesh->number_cells());
  Require(group < d_number_groups);

  return d_factor * d_fixed_source->source(cell, group);
}

//---------------------------------------------------------------------------//
inline double PulsedExternalSource::source(const size_t cell,
                                           const size_t group,
                                           const size_t angle)
{
  // Preconditions
  Require(cell  < d_mesh->number_cells());
  Require(group < d_number_groups);
  Require(angle < d_number_angles);

  return d_factor * d_fixed_source->source(cell, group, angle);
}

//---------------------------------------------------------------------------//
inline void PulsedExternalSource::set_time(const double time)
{

  d_time = time;
  double den = 2.0*std::pow(d_fwhm / 2.354820045030949, 2);
  d_factor = std::exp(-std::pow(d_time-d_peak_time, 2) / den);
  //std::cout << " den=" << den << " factor=" << d_factor << std::endl;
}

} // end namespace detran

#endif // detran_PULSEDEXTERNALSOURCE_I_HH_

//---------------------------------------------------------------------------//
//              end of file PulsedExternalSource.i.hh
//---------------------------------------------------------------------------//
