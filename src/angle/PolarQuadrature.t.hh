//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PolarQuadrature.t.hh
 *  @brief PolarQuadrature member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_POLARQUADRATURE_T_HH_
#define detran_angle_POLARQUADRATURE_T_HH_

#include "angle/PolarQuadrature.hh"
#include "utilities/MathUtilities.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
template <class B>
PolarQuadrature<B>::PolarQuadrature(const size_t number_polar,
                                    bool         normalize)
  : Quadrature(1, 2 * number_polar, B::name())
  , d_number_polar(number_polar)
  , d_cos_theta(number_polar, 0.0)
  , d_sin_theta(number_polar, 0.0)
{
  using detran_utilities::vec_scale;
  using detran_utilities::vec_sum;

  // generate base parameters over [-1, 1]
  B Q;
  Q.build(-1.0, 1.0, 2 * d_number_polar);

  // copy positive half space
  for (size_t i = 0; i < d_number_polar; ++i)
  {
    d_cos_theta[i] = Q.get_x()[i + d_number_polar];
    d_weight[i]    = Q.get_w()[i + d_number_polar];
  }

  // base quadrature stores things in order, but polar quadratures need
  // decreasing polar cosine, so reverse points and weights
//  std::reverse(d_cos_theta.begin(), d_cos_theta.end());
//  std::reverse(d_weight.begin(), d_weight.end());

  // set the sin
  for (size_t i = 0; i < d_number_polar; ++i)
    d_sin_theta[i] = std::sqrt(1.0 - d_cos_theta[i]*d_cos_theta[i]);

  // by default, 1-d problems are along the x-axis
  d_mu = d_cos_theta;

  if (normalize) vec_scale(d_weight, 1.0/vec_sum(d_weight));
}

//----------------------------------------------------------------------------//
template <class B>
typename PolarQuadrature<B>::size_t
PolarQuadrature<B>::number_polar() const
{
  return d_number_polar;
}

//----------------------------------------------------------------------------//
template <class B>
const typename PolarQuadrature<B>::vec_dbl&
PolarQuadrature<B>::sin_theta() const
{
  return d_sin_theta;
}

//----------------------------------------------------------------------------//
template <class B>
const typename PolarQuadrature<B>::vec_dbl&
PolarQuadrature<B>::cos_theta() const
{
  return d_cos_theta;
}

} // end namespace detran_angle

#endif /* detran_angle_POLARQUADRATURE_T_HH_ */

//----------------------------------------------------------------------------//
//              end of PolarQuadrature.t.hh
//----------------------------------------------------------------------------//
