//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ChebyshevU.hh
 *  @brief  ChebyshevU class definition
 *  @author Jeremy Roberts
 *  @date   Jan 23, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_ChebyshevU_HH_
#define detran_orthog_ChebyshevU_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class ChebyshevU
 *  @brief Chebyshev polynomials of the second kind
 */
class ChebyshevU: public ContinuousOrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param   order   Order of the basis
   *   @param   x_0     Lower bound
   *   @param   x_1     Upper bound
   *   @param   x       Grid on which basis is defined
   *   @param   w       Associated weights
   */
  ChebyshevU(const size_t order,
             const vec_dbl &x,
             const vec_dbl &qw,
             const double x_0 = -1.0,
             const double x_1 =  1.0,
             const bool flag = false);

  /// Virtual destructor
  virtual ~ChebyshevU(){}

private :

  double cheby_u(const size_t l, const double x);

};

} // end namespace detran_orthog

#endif // detran_orthog_ChebyshevU_HH_

//---------------------------------------------------------------------------//
//              end of file ChebyshevU.hh
//---------------------------------------------------------------------------//
