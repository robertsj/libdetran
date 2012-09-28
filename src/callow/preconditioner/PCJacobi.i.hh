//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PCJacobi.i.hh
 *  @brief  PCJacobi inline member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

#ifndef PCJACOBI_I_HH_
#define PCJACOBI_I_HH_

namespace callow
{

inline void PCJacobi::apply(Vector &b, Vector &x)
{
  // preconditions
  Require(x.size() == d_P->size());

  // apply x = inv(P)*b
  x.copy(b);
  x.multiply((*d_P));
}


} // end namespace detran

#endif // PCJACOBI_I_HH_ 

//---------------------------------------------------------------------------//
//              end of file PCJacobi.i.hh
//---------------------------------------------------------------------------//
