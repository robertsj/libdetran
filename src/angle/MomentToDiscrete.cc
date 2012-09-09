//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MomentToDiscrete.cc
 * \author Jeremy Roberts
 * \date   Jul 1, 2011
 * \brief  MomentToDiscrete member definitions.
 */
//---------------------------------------------------------------------------//

#include "MomentToDiscrete.hh"

namespace detran_angle
{

// Constructor
MomentToDiscrete::MomentToDiscrete(const size_t legendre_order)
  : d_legendre_order(legendre_order)
  , d_number_moments(1) // \todo save for later use
{
  Require(d_number_moments > 0);
}


} // end namespace slabtran

//---------------------------------------------------------------------------//
//              end of MomentToDiscrete.cc
//---------------------------------------------------------------------------//
