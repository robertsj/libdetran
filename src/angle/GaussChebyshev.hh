//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussChebyshev.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  GaussChebyshev class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_GAUSSCHEBYSHEV_HH_
#define detran_angle_GAUSSCHEBYSHEV_HH_

#include "Quadrature.hh"

namespace detran_angle
{

class GaussChebyshev: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    order       Quadrature order.
   */
  GaussChebyshev(size_t order);

  /// SP constructor
  static SP_quadrature Create(size_t order)
  {
    SP_quadrature p(new GaussChebyshev(order));
    return p;
  }

};

} // end namespace detran_angle


#endif /* detran_angle_GAUSSCHEBYSHEV_HH_ */
