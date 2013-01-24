//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LegendreDTN.cc
 *  @brief  LegendreDTN
 *  @author Jeremy Roberts
 *  @date   Oct 11, 2012
 */
//---------------------------------------------------------------------------//

#include "LegendreDTN.hh"
#include "GenerateGaussLegendre.hh"
#include "GenerateGaussChebyshev.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
LegendreDTN::LegendreDTN(const size_t dim,
                         const size_t na,
                         const size_t np)
  : ProductQuadrature(dim, na, np, "LegendreDTN")
{
  using detran_utilities::pi;

  //-------------------------------------------------------------------------//
  // AZIMUTH QUADRATURE
  //-------------------------------------------------------------------------//

  // temporary arrays
  vec_dbl x(2*na, 0.0);
  vec_dbl w(2*na, 0.0);

  // generate parameters
  generate_gl_parameters(2*na, x, w);

  // fill array
  double w_sum = 0.0;
  for (int i = 0; i < np; ++i)
  {
    d_cos_phi[i] = x[i];
    d_sin_phi[i] = std::sqrt(1.0 - x[i]*x[i]);
    d_azimuth_weight[i] = w[i];
    w_sum += w[i];
  }
  for (int i = 0; i < np; ++i)
  {
    d_azimuth_weight[i] *= 0.5 * pi / w_sum;
  }

  //-------------------------------------------------------------------------//
  // POLAR QUADRATURE
  //-------------------------------------------------------------------------//

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(np, 0.0);
  vec_dbl tmp_wt(np, 0.0);
  generate_gc_parameters(np, tmp_mu, tmp_wt, true);

  w_sum = 0.0;
  for (int i = 0; i < np; i++)
  {
    // Shift the mu to [0, 1] from [-1, 1] and halve the weight
    d_cos_theta[i] = 0.5*tmp_mu[i] + 0.5;
    d_sin_theta[i] = std::sqrt(1.0 - d_cos_theta[i]*d_cos_theta[i]);
    d_polar_weight[i] = tmp_wt[i];
    w_sum += tmp_wt[i];
  }
  for (int i = 0; i < np; ++i)
  {
    d_azimuth_weight[i] *= 1.0 / w_sum;
  }


  //-------------------------------------------------------------------------//
  // PRODUCT QUADRATURE
  //-------------------------------------------------------------------------//

  build_product_quadrature();

}

//---------------------------------------------------------------------------//
LegendreDTN::SP_quadrature
LegendreDTN::Create(const size_t dim,
                     const size_t na,
                     const size_t np)
{
  SP_quadrature p(new LegendreDTN(dim, na, np));
  return p;
}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of file LegendreDTN.cc
//---------------------------------------------------------------------------//
