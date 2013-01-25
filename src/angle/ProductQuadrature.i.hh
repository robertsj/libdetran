//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ProductQuadrature.i.hh
 *  @author robertsj
 *  @date   Oct 10, 2012
 *  @brief  ProductQuadrature.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_PRODUCTQUADRATURE_I_HH_
#define detran_angle_PRODUCTQUADRATURE_I_HH_

namespace detran_angle
{

inline ProductQuadrature::size_t
ProductQuadrature::angle(const size_t a, const size_t p) const
{
  Require(a < d_number_azimuth_octant);
  Require(p < d_number_polar_octant);
  return p + a * d_number_polar_octant;
}

inline ProductQuadrature::size_t
ProductQuadrature::azimuth(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = (angle % d_number_angles_octant) / d_number_polar_octant;
  Ensure(tmp < d_number_azimuth_octant);
  return tmp;
}

inline ProductQuadrature::size_t
ProductQuadrature::polar(const size_t angle) const
{
  Require(angle < d_number_angles_octant);
  size_t tmp = angle % d_number_polar_octant;
  Ensure(tmp < d_number_polar_octant);
  return tmp;
}

inline double ProductQuadrature::sin_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_sin_theta[p];
}

inline double ProductQuadrature::cos_theta(const size_t p) const
{
  Require(p < d_number_polar_octant);
  return d_cos_theta[p];
}

inline double ProductQuadrature::phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_phi[a];
}

inline double ProductQuadrature::sin_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_sin_phi[a];
}

inline double ProductQuadrature::cos_phi(const size_t a) const
{
  Require(a < 2 * d_number_azimuth_octant);
  return d_cos_phi[a];
}

inline ProductQuadrature::size_t
ProductQuadrature::number_azimuths_octant() const
{
  return d_number_azimuth_octant;
}

inline ProductQuadrature::size_t
ProductQuadrature::number_polar_octant() const
{
  return d_number_polar_octant;
}

} // end namespace detran_angle


#endif /* detran_angle_PRODUCTQUADRATURE_I_HH_ */
