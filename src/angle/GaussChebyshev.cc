//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussChebyshev.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  GaussChebyshev class definition.
 */
//---------------------------------------------------------------------------//

#include "GaussChebyshev.hh"
#include "GenerateGaussChebyshev.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
GaussChebyshev::GaussChebyshev(const size_t number_polar_octant,
                               const bool   normalize)
  : Quadrature(1, 2*number_polar_octant, "GaussChebyshev")
{

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(d_number_angles, 0.0);
  vec_dbl tmp_wt(d_number_angles, 0.0);
  generate_gc_parameters(d_number_angles, tmp_mu, tmp_wt, normalize);

  // fill parameters
  for (size_t i = 0; i < d_number_angles_octant; ++i)
  {
    d_mu[i]     = tmp_mu[i];
    d_weight[i] = tmp_wt[i];
  }

}

//---------------------------------------------------------------------------//
GaussChebyshev::SP_quadrature
GaussChebyshev::Create(const size_t number_polar_octant,
                       const bool   normalize)
{
  SP_quadrature p(new GaussChebyshev(number_polar_octant, normalize));
  return p;
}

} // end namespace detran_angle



