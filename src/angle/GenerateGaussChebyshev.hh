//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GenerateGaussChebyshev.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  Gauss-Chebyshev coefficient generator
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_GENERATEGAUSSCHEBYSHEV_HH_
#define detran_angle_GENERATEGAUSSCHEBYSHEV_HH_

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include <cmath>
#include <iostream>

namespace detran_angle
{

/**
 *  @brief Generate Gauss-Chebyshev parameters
 *
 *  The Gauss-Chebyshev quadrature approximates
 *  @f[
 *      \int^{1}_{-1} \frac{f(x)}{\sqrt{1-x^2}} dx
 *        \approx \sum^m_{i=1} W_m f(x_i) \, .
 *  @f]
 *
 *  Now, the integral we want is not actually weighted, and so we
 *  set the true weights to be
 *  @f[
 *     w_i = W_i \sqrt{1-x_i^2} \, .
 *  @f]
 *
 *  The abscissa of the m-point quadrature are the zeros of the
 *  Chebyshev polynomial of degree m + 1.
 *
 *  However, these weights do not sum to two, as we'd expect
 *  when integrating @f$ f(x) = 1 @f$, though the sum does
 *  approach two for high order.  Optionally, the user can
 *  choose to normalize the weights to two.
 *
 *  Reference:
 *    Hildebrand, <i>Introduction to Numerical Analysis</i>
 *
 *  @param m        number of points (i.e. the quadrature order)
 *  @param x        temporary array for abscissa
 *  @param w        temporary array for weights
 *
 */
inline void generate_gc_parameters(detran_utilities::size_t  m,
                                   detran_utilities::vec_dbl &x,
                                   detran_utilities::vec_dbl &w,
                                   bool normalize = false)
{
  Require(x.size() == m);
  Require(w.size() == m);

  const double pi = detran_utilities::pi;

  // Base weight (constant for a particular order)
  double W = pi / (double)m;

  // Generate abscissa
  for (int i = 0; i < m; ++i)
  {
    x[i] = std::cos((2*(i+1)- 1)*pi/(2.0*m));
  }

  // Generate weights
  double w_tot = 0.0;
  for (int i = 0; i < m; ++i)
  {
    w[i] = W * std::sqrt(1.0 - x[i]*x[i]);
    w_tot += w[i];
  }

  // Optionally normalize the weights to 2
  if (normalize)
  {
    for (int i = 0; i < m; ++i)
    {
      w[i] *= 2.0 / w_tot;
    }
  }

}

} // end namespace detran_angle

#endif /* detran_angle_GENERATEGAUSSCHEBYSHEV_HH_ */
