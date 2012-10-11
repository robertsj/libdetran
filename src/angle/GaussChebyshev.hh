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

/**
 *  @class GaussChebyshev
 *  @brief Implements Gauss-Chebyshev quadrature
 *
 *  The Gauss-Chebyshev quadrature approximates
 *  @f[
 *      \int^{1}_{-1} \frac{f(x)}{\sqrt{1-x^2}} dx
 *        \approx \sum^m_{i=1} W_i f(x_i) \, .
 *  @f]
 *  Now, the integral we want is not actually weighted, and so we
 *  set the true weights to be
 *  @f[
 *     w_i = W_i \sqrt{1-x_i^2} \, .
 *  @f]
 *
 *  The abscissa of the m-point quadrature are the zeros of the
 *  Chebyshev polynomial of degree m + 1.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant -- number of abscissa per half space
 *    - gausschebyshev_normalize -- normalize half space weight to 1 [false]
 */
class GaussChebyshev: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_polar_octant    number of angles per half space
   *  @param normalize              normalize half space weights to 1
   */
  GaussChebyshev(const size_t number_polar_octant,
                 const bool   normalize = false);

  /// SP constructor
  static SP_quadrature Create(const size_t number_polar_octant,
                              const bool   normalize = false);

};

} // end namespace detran_angle


#endif /* detran_angle_GAUSSCHEBYSHEV_HH_ */
