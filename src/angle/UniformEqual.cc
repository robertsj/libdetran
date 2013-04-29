//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UniformEqual.cc
 * \author robertsj
 * \date   May 22, 2012
 * \brief  UniformEqual member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#include "UniformEqual.hh"
#include <iostream>
#include <cmath>
#include <cstdio>

namespace detran_angle
{

UniformEqual::UniformEqual(size_t order, size_t dim)
  : Quadrature(dim,
               std::pow((float)2, (int)dim) * order * order,
               "UniformEqual")
{

  Require(order >= 1);

  const double pi = detran_utilities::pi;

  // Total number of polar angles
  size_t num_xi = 2 * order;

  // Total number of azimuthal angles
  size_t num_omega = 2 * num_xi;

  Ensure(num_xi * num_omega == d_number_angles);
  Ensure(num_omega % 4 == 0);

  // Common values below
  size_t num_omega_4   = num_omega / 4; // omega per octant
  size_t num_xi_2      = num_xi / 2;    // xi per octant
  double two_num_xi = 2.0 / num_xi;  //

  // temporary arrays to hold the positive octant values
  vec_dbl tmp_xi(num_xi_2 + 1, 0.0); // +1 to handle order=2 case
  // uniform xi density
  if (1)
  {
    // uniform values of the cosine of the polar angle
    for (size_t i = 0; i < num_xi_2 + 1; i++) // +1 to handle order=2 case
      tmp_xi[i] = two_num_xi * (i + 0.5 + 0.5 * num_xi) - 1.0;
  }
  else
  {
    THROW("Not yet implemented");
  }

  // width of polar angle cosine mesh
  double del_xi = tmp_xi[1] - tmp_xi[0];

  double coef = 4.0 * pi / (d_number_angles_octant * 8 * del_xi);
  double tmp = 0;
  double cos_omega = 0;
  double sin_omega = 0;

  // Uniform weight
  double wt = detran_utilities::four_pi / d_number_angles;

  // compute mu and eta
  size_t k = 0;
  for (size_t j = 0; j < num_omega_4; j++)
  {
    tmp = coef * (j + 0.5);
    cos_omega = std::cos(tmp);
    sin_omega = std::sin(tmp);
    for (size_t i = 0; i < num_xi_2; i++)
    {
      d_xi[k]  = tmp_xi[i];
      d_mu[k]  = cos_omega * std::sqrt(1.0 - tmp_xi[i] * tmp_xi[i]);
      d_eta[k] = sin_omega * std::sqrt(1.0 - tmp_xi[i] * tmp_xi[i]);
      d_weight[k] = wt;
      k++;
    }
  }
  Ensure(k == d_number_angles_octant);

}

} // end namespace detran_angle

//---------------------------------------------------------------------------//
//              end of UniformEqual.cc
//---------------------------------------------------------------------------//


