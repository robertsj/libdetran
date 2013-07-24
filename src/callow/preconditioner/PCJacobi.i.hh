//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCJacobi.i.hh
 *  @brief PCJacobi inline member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PCJACOBI_I_HH_
#define callow_PCJACOBI_I_HH_

namespace callow
{

inline void PCJacobi::apply(Vector &b, Vector &x)
{
  Require(x.size() == d_P->size());

  // apply x = inv(P)*b
  x.copy(b);
  x.multiply((*d_P));
}

} // end namespace detran

#endif // callow_PCJACOBI_I_HH_

//----------------------------------------------------------------------------//
//              end of file PCJacobi.i.hh
//----------------------------------------------------------------------------//
