//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ChebyshevU.hh
 *  @brief ChebyshevU class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_ChebyshevU_HH_
#define detran_orthog_ChebyshevU_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class ChebyshevU
 *  @brief Chebyshev polynomials of the second kind
 */
class ORTHOG_EXPORT ChebyshevU: public ContinuousOrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Constructor
  ChebyshevU(const Parameters &p);

  /// Virtual destructor
  virtual ~ChebyshevU(){}

private :

  double cheby_u(const size_t l, const double x);

};

} // end namespace detran_orthog

#endif // detran_orthog_ChebyshevU_HH_

//----------------------------------------------------------------------------//
//              end of file ChebyshevU.hh
//----------------------------------------------------------------------------//
