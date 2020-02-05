//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  AzimuthalQuadrature.t.hh
 *  @brief AzimuthalQuadrature member definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_AZIMUTHALQUADRATURE_T_HH_
#define detran_angle_AZIMUTHALQUADRATURE_T_HH_

#include "angle/AzimuthalQuadrature.hh"
#include "utilities/Constants.hh"
#include "utilities/MathUtilities.hh"

namespace detran_angle
{

//----------------------------------------------------------------------------//
template <class B>
AzimuthalQuadrature<B>::AzimuthalQuadrature(const size_t number_azimuth,
                                            bool         normalize)
  : d_number_azimuth(number_azimuth)
  , d_phi(number_azimuth, 0.0)
  , d_cos_phi(number_azimuth, 0.0)
  , d_sin_phi(number_azimuth, 0.0)
  , d_weight(number_azimuth, 0.0)
{
  using detran_utilities::vec_scale;
  using detran_utilities::vec_sum;

  // generate base parameters over [0, pi/2]
  B Q;
  Q.build(0.0, std::acos(0.0), d_number_azimuth);

  // copy from base
  d_phi = Q.get_x();
  for (size_t a = 0; a < d_number_azimuth; ++a)
  {
    d_cos_phi[a] = std::cos(d_phi[a]);
    d_sin_phi[a] = std::sin(d_phi[a]);
  }
  d_weight = Q.get_w();

  // note
  if (normalize) vec_scale(d_weight, std::acos(0.0)/vec_sum(d_weight));
}

//----------------------------------------------------------------------------//
template <class B>
typename AzimuthalQuadrature<B>::size_t
AzimuthalQuadrature<B>::number_azimuth() const
{
  return d_number_azimuth;
}

//----------------------------------------------------------------------------//
template <class B>
const typename AzimuthalQuadrature<B>::vec_dbl&
AzimuthalQuadrature<B>::phi() const
{
  return d_phi;
}

//----------------------------------------------------------------------------//
template <class B>
const typename AzimuthalQuadrature<B>::vec_dbl&
AzimuthalQuadrature<B>::cos_phi() const
{
  return d_cos_phi;
}

//----------------------------------------------------------------------------//
template <class B>
const typename AzimuthalQuadrature<B>::vec_dbl&
AzimuthalQuadrature<B>::sin_phi() const
{
  return d_sin_phi;
}

//----------------------------------------------------------------------------//
template <class B>
const typename AzimuthalQuadrature<B>::vec_dbl&
AzimuthalQuadrature<B>::weights() const
{
  return d_weight;
}

} // end namespace detran_angle

#endif /* detran_angle_AZIMUTHALQUADRATURE_T_HH_ */

//----------------------------------------------------------------------------//
//              end of file AzimuthalQuadrature.t.hh
//----------------------------------------------------------------------------//
