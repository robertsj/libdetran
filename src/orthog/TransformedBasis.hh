//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  TransformedBasis.hh
 *  @brief TransformedBasis class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_TRANSFORMEDBASIS_HH_
#define detran_orthog_TRANSFORMEDBASIS_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class TransformedBasis
 *  @brief Discrete basis based on user-defined transformation of zeroth moment
 */
class ORTHOG_EXPORT TransformedBasis: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor.  The user-defined zeroth order basis lives in p.x.
  TransformedBasis(const Parameters &p);

  /// Virtual destructor
  virtual ~TransformedBasis(){}

  //AUTO_REGISTER(TransformedBasis, OrthogonalBasis)
};

} // end namespace detran_orthog

#endif // detran_orthog_TRANSFORMEDBASIS_HH_

//----------------------------------------------------------------------------//
//              end of file TransformedBasis.hh
//----------------------------------------------------------------------------//
