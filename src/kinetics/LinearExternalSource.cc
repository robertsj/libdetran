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
                                           vec_source     sources)
  : Base(number_groups, mesh, SP_quadrature(0))
  , d_times(times)
  , d_sources(sources)
{
   // Preconditions
   Require(d_times.size() > 0);
   Require(d_times.size() == d_sources.size());

   d_number_times = d_times.size();
}

//---------------------------------------------------------------------------//
LinearExternalSource::SP_externalsource
LinearExternalSource::Create(const size_t   number_groups,
                             SP_mesh        mesh,
                             vec_dbl        times,
                             vec_source     sources)
{
  SP_externalsource p(new LinearExternalSource(number_groups, mesh,
                                               times, sources));
  return p;
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file LinearExternalSource.cc
//---------------------------------------------------------------------------//
