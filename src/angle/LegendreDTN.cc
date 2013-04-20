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
  vec_dbl x(na, 0.0);
  vec_dbl w(na, 0.0);

  // generate parameters
  generate_gl_parameters(na, x, w);

  for (int i = 0; i < na; ++i)
  {
    std::cout << " i = " << i << " x = " << x[i] << std::endl;
  }

  // fill array
  double w_sum = 0.0;
  for (int i = 0; i < na; ++i)
  {
    size_t j     = na - i - 1;
    d_phi[j]     = 0.25*pi*(x[i] + 1.0);
    d_cos_phi[j] = std::cos(d_phi[j]);
    d_sin_phi[j] = std::sin(d_phi[j]);
    d_azimuth_weight[j] = w[i];
    w_sum += w[i];
  }
  for (int i = 0; i < na; ++i)
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
    // Put in order.
    size_t j = np - i - 1;
    // Shift the mu to [0, 1] from [-1, 1] and halve the weight
    d_cos_theta[j] = 0.5*tmp_mu[i] + 0.5;
    d_sin_theta[j] = std::sqrt(1.0 - d_cos_theta[j]*d_cos_theta[j]);
    d_polar_weight[j] = tmp_wt[i];
    w_sum += tmp_wt[i];
  }
  for (int i = 0; i < np; ++i)
  {
    d_polar_weight[i] *= 1.0 / w_sum;
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
