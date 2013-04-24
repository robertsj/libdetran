//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CLP.hh
 *  @brief  CLP class definition
 *  @author Jeremy Roberts
 *  @date   Jan 8, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_CLP_HH_
#define detran_orthog_CLP_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class CLP
 *  @brief Continuous Legendre polynomial basis
 */
class ORTHOG_EXPORT CLP: public ContinuousOrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  CLP(const size_t order,
      const vec_dbl &x,
      const vec_dbl &qw,
      const double x_0 = -1.0,
      const double x_1 = 1.0);

  /// Virtual destructor
  virtual ~CLP(){}

};


} // end namespace detran_orthog

#endif // detran_orthog_CLP_HH_

//---------------------------------------------------------------------------//
//              end of file CLP.hh
//---------------------------------------------------------------------------//
