//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DPN.cc
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  DPN class definition.
 */
//---------------------------------------------------------------------------//

#include "DPN.hh"
#include "GenerateGaussLegendre.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
DPN::DPN(const size_t number_polar_octant)
  : Quadrature(1, 2*number_polar_octant, "DPN")
{

  // generate the parameters on-the-fly.
  vec_dbl tmp_mu(d_number_angles_octant, 0.0);
  vec_dbl tmp_wt(d_number_angles_octant, 0.0);
  generate_gl_parameters(d_number_angles_octant, tmp_mu, tmp_wt);

  for (int i = 0; i < d_number_angles_octant; i++)
  {
    // Shift the mu to [0, 1] from [-1, 1]
    d_mu[i]     = 0.5*tmp_mu[i] + 0.5;
    d_weight[i] = 0.5*tmp_wt[i];
  }

}

//---------------------------------------------------------------------------//
DPN::SP_quadrature
DPN::Create(const size_t number_polar_octant)
{
  SP_quadrature p(new DPN(number_polar_octant));
  return p;
}

} // end namespace detran_angle


