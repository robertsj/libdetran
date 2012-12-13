//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MomentToDiscrete.cc
 *  @author Jeremy Roberts
 *  @date   Jul 1, 2011
 *  @brief  MomentToDiscrete member definitions.
 */
//---------------------------------------------------------------------------//

#include "MomentToDiscrete.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
MomentToDiscrete::MomentToDiscrete(SP_momentindexer indexer)
  : d_indexer(indexer)
{
  // Preconditioner
  Require(d_indexer);

  d_legendre_order = d_indexer->legendre_order();
  d_number_moments = d_indexer->number_moments();
}


} // end namespace slabtran

//---------------------------------------------------------------------------//
//              end of MomentToDiscrete.cc
//---------------------------------------------------------------------------//
