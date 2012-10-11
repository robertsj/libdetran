//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DPN.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  DPN class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_DPN_HH_
#define detran_angle_DPN_HH_

// Angle headers
#include "Quadrature.hh"

namespace detran_angle
{

/**
 *  @class DPN
 *  @brief Implements the double GaussLegendre quadrature
 *
 *  This is identical to Gauss-Legendre quadrature except that the
 *  quadrature parameters are based on the integral
 *  @f[
 *      \int^{1}_{0} f(x) dx
 *        \approx \sum^m_{i=1} W_m f(x_i) \, .
 *  @f]
 *  Doing this effectively treats the left and right half spaces
 *  independently, which can work well for functions that have completely
 *  different behavior in each half space but within those half spaces, the
 *  behavior is polynomial-like.
 *
 *  The resulting abscissa are those of the GL quadrature but scaled
 *  and shifted.
 *
 */
class DPN: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param order       Quadrature order (total number of polar angles)
   */
  DPN(size_t order);

  /// SP constructor
  static SP_quadrature Create(size_t order)
  {
    SP_quadrature p(new DPN(order));
    return p;
  }

};

} // end namespace detran_angle

#endif /* detran_angle_DPN_HH_ */
