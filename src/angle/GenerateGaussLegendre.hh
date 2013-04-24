//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GenerateGaussLegendre.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  Gauss-Legendre coefficient generator
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_GENERATEGAUSSLEGENDRE_HH_
#define detran_angle_GENERATEGAUSSLEGENDRE_HH_

#include "utilities/Definitions.hh"
#include <cmath>

namespace detran_angle
{

/**
 *  @brief Generate Gauss-Legendre parameters.
 *
 *  The Gauss-Legendre quadrature approximates
 *  @f[
 *      \int^{1}_{-1} f(x) dx
 *        \approx \sum^m_{i=1} W_m f(x_i) \, .
 *  @f]
 *
 *  The abscissa @f$ x_i @f$ of the m-point quadrature
 *  turn out to be the zeroes of the
 *  Legendre polynomial of degree @f$ m + 1 @f$.  The weights
 *  can be found numerically, as done in this implementation,
 *  which is a modified version of the function
 *  gauleg from <B> Numerical Recipes in C</B>
 *  by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
 *  & B.P. Flannery, (Cambridge Univ. Press)
 *
 *  @param m        number of points (i.e. the quadrature order)
 *  @param x        temporary array for abscissa
 *  @param w        temporary array for weights
 *
 */
inline void generate_gl_parameters(detran_utilities::size_t  m,
                                   detran_utilities::vec_dbl &x,
                                   detran_utilities::vec_dbl &w)
{
  Require(x.size() == m);
  Require(w.size() == m);

  size_t j, i;
  double z1, z, pp, p3, p2, p1;

  // The roots are symmetric, so we only find half of them.
  for (i = 1; i <= m; i++)
  { /* Loop over the desired roots. */
    z = std::cos(detran_utilities::pi * (i - 0.25) / (m + 0.5));
    // Starting with the above approximation to the ith root, we enter
    // the main loop of refinement by Newon's method.
    int count = 0;
    do
    {
      p1 = 1.0;
      p2 = 0.0;
      // Recurrence to get Legendre polynomial.
      for (j = 1; j <= m; j++)
      {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
      }
      // p1 is now the desired Legendre polynomial. We next compute
      // pp, its derivative, by a standard relation involving also
      // p2, the polynomial of one lower order.
      pp = m * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp; // <-- Newon's method
    } while (std::abs(z - z1) > 1e-15 && ++count < 200);
    // Store the root and compute the weight.
    x[i - 1] = z;
    w[i - 1] = 2.0 / ((1.0 - z * z) * pp * pp);
  }
}

} // end namespace detran_angle

#endif /* detran_angle_GENERATEGAUSSLEGENDRE_HH_ */
