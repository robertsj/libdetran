//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   UniformEqual.cc
 * \author robertsj
 * \date   May 22, 2012
 * \brief  UniformEqual member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "UniformEqual.hh"

// System
#include <iostream>
#include <cmath>
#include <cstdio>

namespace detran
{

UniformEqual::UniformEqual(int order, int dim)
  : Quadrature(order,
               dim,
               std::pow(2, dim) * order * order,
               "UniformEqual")
{

  Require(order >= 1);

  // Total number of polar angles
  int num_xi = 2 * order;

  // Total number of azimuthal angles
  int num_omega = 2 * num_xi;

  Ensure(num_xi * num_omega == d_number_angles);
  Ensure(num_omega % 4 == 0);

  // Common values below
  int num_omega_4   = num_omega / 4; // omega per octant
  int num_xi_2      = num_xi / 2;    // xi per octant
  double two_num_xi = 2.0 / num_xi;  //

  // temporary arrays to hold the positive octant values
  vec_dbl tmp_xi(num_xi_2 + 1, 0.0); // +1 to handle order=2 case
  // uniform xi density
  if (1)
  {
    // uniform values of the cosine of the polar angle
    for (int i = 0; i < num_xi_2 + 1; i++) // +1 to handle order=2 case
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
  double wt = four_pi / d_number_angles;

  // compute mu and eta
  int k = 0;
  for (int j = 0; j < num_omega_4; j++)
  {
    tmp = coef * (j + 0.5);
    cos_omega = std::cos(tmp);
    sin_omega = std::sin(tmp);
    for (int i = 0; i < num_xi_2; i++)
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

void UniformEqual::display() const
{

    using std::cout;
    using std::endl;
    using std::printf;

    cout << endl;
    cout << d_name << " abscissa and weights: " << endl << endl;
    cout << "   m            mu                 eta                xi                 wt       " << endl;
    cout << "  ---   -----------------  -----------------  -----------------  -----------------" << endl;
    double weight_sum = 0.0;
    for ( int ix = 0; ix < d_number_angles_octant; ++ix )
    {
        printf ("%4i    %16.13f   %16.13f   %16.13f   %16.13f   \n",
                ix, d_mu[ix], d_eta[ix], d_xi[ix], d_weight[ix] );
        weight_sum += d_weight[ix];
    }
    cout << endl << "  The sum of the weights is " << weight_sum << endl;
    cout << endl;

}

} // end namespace detran

