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
 *  Relevant database parameters:
 *    - quad_number_polar_octant -- number of abscissa per half space
 *    - gausschebyshev_normalize -- normalize half space weight to 1 [false]
 */
class ANGLE_EXPORT DTN: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_polar_octant    number of angles per half space
   *  @param normalize              normalize half space weights to 1
   */
  DTN(const size_t number_polar_octant,
      const bool   normalize = false);

  /// SP constructor
  static SP_quadrature Create(const size_t number_polar_octant,
                              const bool   normalize = false);

};

} // end namespace detran_angle

#endif /* detran_angle_DTN_HH_ */
