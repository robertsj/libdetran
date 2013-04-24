//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LinearExternalSource.cc
 *  @brief  LinearExternalSource
 *  @author Jeremy Roberts
 *  @date   Nov 16, 2012
 */
//---------------------------------------------------------------------------//

#include "LinearExternalSource.hh"

namespace detran
{

//---------------------------------------------------------------------------//
LinearExternalSource::LinearExternalSource(const size_t   number_groups,
                                           SP_mesh        mesh,
                                           vec_dbl        times,
                                           vec_source     sources,
                                           bool           discrete)
  : Base(number_groups, mesh, SP_quadrature(0), discrete)
  , d_times(times)
  , d_sources(sources)
{
   // Preconditions
   Require(d_times.size() > 0);
   Require(d_times.size() == d_sources.size());

   d_number_times = d_times.size();
}

//---------------------------------------------------------------------------//
LinearExternalSource::SP_tdsource
LinearExternalSource::Create(const size_t   number_groups,
                             SP_mesh        mesh,
                             vec_dbl        times,
                             vec_source     sources,
                             bool           discrete)
{
  SP_tdsource p(new LinearExternalSource(number_groups, mesh,
                                         times, sources, discrete));
  return p;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file LinearExternalSource.cc
//---------------------------------------------------------------------------//
