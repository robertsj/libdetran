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
 *  Relevant database parameters:
 *    - quad_number_polar_octant -- number of abscissa per half space
 */
class ANGLE_EXPORT DPN: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_polar_octant  Number of polar angles per octant
   */
  DPN(const size_t number_polar_octant);

  /// SP constructor
  static SP_quadrature Create(const size_t number_polar_octant);


};

} // end namespace detran_angle

#endif /* detran_angle_DPN_HH_ */
