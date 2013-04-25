//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Quadrature.i.hh
 *  @brief  Quadrature inline member definitions
 *  @author Jeremy Roberts
 *  @date   Sep 5, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_QUADRATURE_I_HH_
#define detran_angle_QUADRATURE_I_HH_

namespace detran_angle
{

inline Quadrature::size_t Quadrature::number_angles() const
{
  return d_number_angles;
}

inline Quadrature::size_t Quadrature::number_octants() const
{
  return d_number_octants;
}

inline Quadrature::size_t Quadrature::number_angles_octant() const
{
  return d_number_angles_octant;
}

//inline Quadrature::size_t Quadrature::number_azimuths_octant() const
//{
//  THROW("A PRODUCT QUADRATURE MUST IMPLEMENT THIS");
//  return 0;
//}
//
//inline Quadrature::size_t Quadrature::number_polar_octant() const
//{
//  THROW("A PRODUCT QUADRATURE MUST IMPLEMENT THIS");
//  return 0;
//}

inline Quadrature::size_t Quadrature::index(const size_t o, const size_t a)
{
  Require(o < d_number_octants);
  Require(a < d_number_angles_octant);
  size_t angle = a + o * d_number_angles_octant;
  Ensure(angle < d_number_angles);
  return angle;
}

//inline Quadrature::size_t Quadrature::angle(const size_t a, const size_t p) const
//{
//  THROW("NOT A PRODUCT QUADRATURE");
//  return 0;
//}

inline const Quadrature::vec_dbl& Quadrature::weights() const
{
  return d_weight;
}

inline const Quadrature::vec_dbl& Quadrature::cosines(const size_t dir) const
{
  if (dir == MU)
  {
    return d_mu;
  }
  else if (dir == ETA)
  {
    Require(d_dimension > 1);
    return d_eta;
  }
  else
  {
    Require(d_dimension > 1);
    return d_xi;
  }
}

inline double Quadrature::weight(const size_t a) const
{
  Require(a < d_number_angles_octant);
  return d_weight[a];
}

inline double Quadrature::mu(const size_t o, const size_t a) const
{
  Require(a < d_number_angles_octant);
  Require(o < d_number_octants);
  return d_octant_sign[o][MU] * d_mu[a];
}

inline double Quadrature::eta(const size_t o, const size_t a) const
{
  Require(a < d_number_angles_octant);
  Require(o < d_number_octants);
  return d_octant_sign[o][ETA] * d_eta[a];
}

inline double Quadrature::xi(const size_t o, const size_t a) const
{
  Require(a < d_number_angles_octant);
  Require(o < d_number_octants);
  return d_octant_sign[o][XI] * d_xi[a];
}

inline double Quadrature::angular_norm(const size_t d)
{
  Require(d < 4);
  if (d == 1)
    return 0.5;
  else
    return detran_utilities::inv_four_pi;
}

inline bool Quadrature::valid_index(const size_t o, const size_t a) const
{
  if (o < d_number_octants && a < d_number_angles_octant)
  {
    return true;
  }
  else
  {
    return false;
  }
}

} // end namespace detran_angle

#endif // detran_angle_QUADRATURE_I_HH_

//---------------------------------------------------------------------------//
//              end of file Quadrature.i.hh
//---------------------------------------------------------------------------//
