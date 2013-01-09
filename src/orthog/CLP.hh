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
class CLP: public ContinuousOrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param   order   Order of the basis
   *   @param   x       Grid on which basis is defined
   *   @param   dx      Cell widths
   */
  CLP(const size_t order,
      const vec_dbl &x,
      const vec_dbl &dx);

  /// Virtual destructor
  virtual ~CLP(){}

};


} // end namespace detran_orthog

#endif // detran_orthog_CLP_HH_

//---------------------------------------------------------------------------//
//              end of file CLP.hh
//---------------------------------------------------------------------------//
