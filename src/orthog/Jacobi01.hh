//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Jacobi01.hh
 *  @brief  Jacobi01 class definition
 *  @author Jeremy Roberts
 *  @date   Jan 21, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_Jacobi01_HH_
#define detran_orthog_Jacobi01_HH_

#include "ContinuousOrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class Jacobi01
 *  @brief Continuous Jacobi_01 polynomial basis
 *
 *  Implements the Jacobi polynomial with alpha = 0, beta = 1.
 */
class Jacobi01: public ContinuousOrthogonalBasis
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
  Jacobi01(const size_t order,
           const vec_dbl &x,
           const vec_dbl &qw,
           const double x_0 = -1.0,
           const double x_1 =  1.0);

  /// Virtual destructor
  virtual ~Jacobi01(){}

};


} // end namespace detran_orthog

#endif // detran_orthog_Jacobi01_HH_

//---------------------------------------------------------------------------//
//              end of file Jacobi01.hh
//---------------------------------------------------------------------------//
