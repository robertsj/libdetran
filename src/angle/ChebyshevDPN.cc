//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ChebyshevDPN.cc
 *  @brief ChebyshevDPN member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "ChebyshevDPN.hh"
#include "GenerateGaussLegendre.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
ChebyshevDPN::ChebyshevDPN(const size_t dim, const size_t na, const size_t np)
  : ProductQuadrature(dim, na, np, "ChebyshevDPN")
{

  //--------------------------------------------------------------------------//
  // AZIMUTH QUADRATURE
  //--------------------------------------------------------------------------//

  // Chebyshev points are just equally spaced phi's in the first quadrant
  using detran_utilities::pi;
  for (size_t i = 0; i < na; ++i)
  {
    d_phi[i] = 0.25 * (2*i + 1)*pi / na;
    d_cos_phi[i] = std::cos(d_phi[i]);
    d_sin_phi[i] = std::sin(d_phi[i]);
    d_azimuth_weight[i] = 0.5 * pi / na;
  }

  //--------------------------------------------------------------------------//
  // POLAR QUADRATURE
  //--------------------------------------------------------------------------//

  // temporary arrays
  vec_dbl x(np, 0.0);
  vec_dbl w(np, 0.0);

  // generate parameters
  generate_gl_parameters(np, x, w);

  // fill array
  for (size_t i = 0; i < np; ++i)
  {
    // Put in order.
    size_t j = np - i - 1;
    d_cos_theta[j]    = 0.5*x[i] + 0.5; // shift and scale
    d_sin_theta[j]    = std::sqrt(1.0 - d_cos_theta[j]*d_cos_theta[j]);
    d_polar_weight[j] = 0.5*w[i];
  }

  //--------------------------------------------------------------------------//
  // PRODUCT QUADRATURE
  //--------------------------------------------------------------------------//

  build_product_quadrature();

}

//----------------------------------------------------------------------------//
ChebyshevDPN::SP_quadrature
ChebyshevDPN::Create(const size_t dim, const size_t na, const size_t np)
{
  SP_quadrature p(new ChebyshevDPN(dim, na, np));
  return p;
}

} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of file ChebyshevDPN.cc
//----------------------------------------------------------------------------//
