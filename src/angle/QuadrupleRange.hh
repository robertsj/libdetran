//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   QuadrupleRange.hh
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 *  @brief  QuadrupleRange class definition.
 */
//---------------------------------------------------------------------------//

#ifndef QUADRUPLERANGE_HH_
#define QUADRUPLERANGE_HH_

// Detran
#include "Quadrature.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
/**
 * @class QuadrupleRange
 * @brief Quadruple Range quadrature class.
 *
 * For two-dimensional transport, Abu-Shumays developed several quadrature
 * sets that are in some sense analogs of the DPn quadrature.  The sets he
 * developed are product quadratures, which means that polar and azimuthal
 * angles are independent. To illustrate, consider a function
 * \f$ f(\Omega_x,\Omega_y) \f$. Using a product quadrature, the integral
 * is approximated as
 * \f[
 *    \int_{\Omega} f(\Omega_x,\Omega_y) d\Omega \approx
 *    \sum^{4N_{\phi}}_{i=1} \sum^{N_{\theta}}_{j=1}
 *        w_i w_j f(\sin \theta_j \cos \phi_i, \sin \theta_j \sin \phi_i ) \, ,
 * \f]
 * If a function \f$ f \f$ is a polynomial in \f$\Omega_x \f$
 * and \f$\Omega_y \f$ of degree \f$ L \f$, a quadrature set that can
 * exactly integrate \f$ f \f$ is denoted order \f$ L \f$.
 *
 * Here, we implement the quadruple range quadrature, which Abu-Shumays
 * estimated would be well-suited for ``highly heterogeneous problems with
 * numerous corner singularities''.   The quadrature essentially provides a
 * good approximation for angular fluxes that are represented by distinct
 * polynomials in \f$\Omega_x \f$ and \f$\Omega_y \f$ within each quadrant.
 * We expect such a quadrature to be well-suited for response function
 * generation, where highly anisotropic fluxes will result due to incident
 * flux conditions restricted to single surfaces of a region.
 *
 * The polar angle quadrature is Abu-Shumays' double-range quadrature, labeled
 * \f$ K^d_L \f$, which satisfies
 * \f[
         \int^{\pi}_{0} d\theta \sin \theta \sin^{j} \! \theta =
          \int^{1}_{0} x^j \frac{2xdx}{\sqrt{1-x^2}} \, ,
 * \f]
 * where \f$ x = \sin \theta \f$. To represent terms for \f$ j = 0 \ldots L \f$
 * requires \f$ L+1 \f$ equations, sufficient for determining \f$ (L+1)/2 \f$
 * points and \f$ (L+1)/2 \f$ weights.  Note, this quadrature is nothing more
 * that a
 * Gaussian quadrature with weight \f$ w(x)=2x/\sqrt{1-x^2} \f$ over the range
 * \f$ 0 \leq x \leq 1 \f$.  Note furthermore that \f$ L = 2N_{\theta}-1 \f$.
 *
 * The azimuthal angle quadrature is Abu-Shumays' symmetric quadruple-range
 * quadrature, labeled \f$ I^{q0}_{L} \f$.  This quadrature satisfies
 * \f[
        \int^{\pi/2}_{0} d\phi \cos^n \! \phi \sin^m \! \phi
         \approx \sum^{N_{\phi}}_{i=1} w_i \cos^n \! \phi_i \sin^m \! \phi_i \, .
 * \f]
 * subject to the symmetry constraints
 * \f[
        \phi_{N_{\phi}+1-i} = \frac{\pi}{2} - \phi_i \, , \, \, \, \, \,
         i = 1 \ldots \frac{N_{\phi}+1}{2} \,
 * \f]
 * and
 * \f[
       w_{N_{\phi}+1-i} = w_i \, , \, \, \, \, \, i = 1 \ldots \frac{N_{\phi}+1}{2} \, .
 * \f]
 * Given these constraints, the moment equations are not all independent, and it
 * suffices to consider only integrands of form
 * \f[
       \sin^{4k} \!  \phi\, , \,\, \sin^{2k+1} \! \phi \, , \, \,
        \sin^{4k+1} \! \phi \cos  \phi \, , \,\,\,\,\, k = 0,1\ldots \, .
 * \f]
 *
 * Suppose we select the number of points \f$ N_{\phi}\f$, which requires
 * \f$ 2N_{\phi} \f$ equations.  The symmetry constraints account for half of these equations.
 * The remaining \f$ N_{\phi} \f$ degrees of freedom allow for exact integration of
 * polynomials of degree \f$ L = N_{\phi}-1 \f$.  Thus, for a given \f$ L \f$, we require twice
 * as many azimuthal points as polar points.  When each set has the same \f$ L \f$, they
 * are said to be compatible, i.e. they have the same order of accuracy.  Moreover,
 * it gives a decisive way to choose the number of polar angles for a given choice
 * of the azimuthal angle count.
 *
 * For the quadratures as implemented, we use double-precision values as
 * generated in Maple and Mathematica in tables from which the full set
 * is generated.  In the future, we could also hard code the full
 * set (via Maple's codegen feature).
 *
 * The constructor uses the number of azimuthal angles, represented
 * by the sn_order parameter.
 *
 * References:
 *
 * Abu-Shumays, I.K. <em>Nuclear Science and Engineering</em>
 *     <b>64</b>, 299-316 (1977).
 *
 *
 *
 */
/**
 *  \example angle/test/test_QuadrupleRange.cc
 *
 *  Test of class QuadrupleRange
 */
//---------------------------------------------------------------------------//

class QuadrupleRange : public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    order       Quadrature order.
   */
  QuadrupleRange(const size_t order, const size_t dim = 2);

  /// SP constructor
  static SP_quadrature Create(const size_t order, const size_t dim = 2);

private:

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /**
   *  @brief Return cosine of \f$ \theta \f$.
   *  @param    N   number of azimuthal angles per quadrant
   *  @param    i   0 for angle and 1 for weight
   *  @param    j   index of requested point
   */
  double get_theta(size_t N, size_t i, size_t j);

  /**
   *  @brief Return sine of \f$ \phi \f$.
   *  @param    N   number of polar angles per quadrant
   *  @param    i   0 for angle and 1 for weight
   *  @param    j   index of requested point
   */

  double get_phi(size_t N, size_t i, size_t j);

};

} // end namespace detran_angle

#endif /* QUADRUPLERANGE_HH_ */

//---------------------------------------------------------------------------//
//              end of QuadrupleRange.hh
//---------------------------------------------------------------------------//
