//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DTN.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  DTN class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_DTN_HH_
#define detran_angle_DTN_HH_

// Angle headers
#include "Quadrature.hh"

namespace detran_angle
{

/**
 *  @class DTN
 *  @brief Implements the double GaussChebyshev quadrature
 *
 *  This is identical to Gauss-Chebyshev quadrature except that the
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
 *  The resulting abscissa are those of the GC quadrature but scaled
 *  and shifted.
 *
 */
class DTN: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param order       Quadrature order (total number of polar angles)
   */
  DTN(size_t order);

  /// SP constructor
  static SP_quadrature Create(size_t order)
  {
    SP_quadrature p(new DTN(order));
    return p;
  }

};

} // end namespace detran_angle

#endif /* detran_angle_DTN_HH_ */
