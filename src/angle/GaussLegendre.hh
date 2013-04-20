//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussLegendre.hh
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 *  @brief  GaussLegendre class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_GAUSSLEGENDRE_HH_
#define detran_angle_GAUSSLEGENDRE_HH_

#include "Quadrature.hh"

namespace detran_angle
{

/**
 *  @class GaussLegendre
 *  @brief Implements Gauss-Legendre quadrature
 *
 *  Gauss-Legendre quadrature approximates
 *  @f[
 *      \int^{1}_{-1} f(x) dx
 *        \approx \sum^m_{i=1} w_i f(x_i) \, .
 *  @f]
 *  The abscissa @f$ x_i @f$ are the zeros of the Legendre
 *  polynomial of degree @f$ m + 1 @f$.  The weights are
 *  then computed numerically.  Notably, an m-point G-L
 *  quadrature exactly integrates polynomials of degree
 *  less than or equal to 2m - 1.  In general, this is
 *  the best that can be achieved when integrating
 *  polynomial functions.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant -- number of abscissa per half space
 */
class GaussLegendre: public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param number_polar_octant  Number of polar angles per octant
   */
  GaussLegendre(const size_t number_polar_octant);

  /// SP constructor
  static SP_quadrature Create(const size_t number_polar_octant);

};

} // end namespace detran_angle

#endif /* detran_angle_GAUSSLEGENDRE_HH_ */

//---------------------------------------------------------------------------//
//              end of GaussLegendre.hh
//---------------------------------------------------------------------------//
