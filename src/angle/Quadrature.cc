//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Quadrature.cc
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "Quadrature.hh"

// System 
#include <iostream>
#include <cstdio>
#include <cmath>

namespace detran
{

// Constructor
Quadrature::Quadrature(int order,
                       int dim,
                       int number_angles,
                       std::string name)
  : d_dimension(dim)
  , d_order(order)
  , d_number_angles(number_angles)
  , d_number_octants(std::pow(2, dim))
  , d_number_angles_octant(number_angles/d_number_octants)
  , d_weight(d_number_angles_octant, 0.0)
  , d_mu(d_number_angles_octant, 0.0)
  , d_eta(d_number_angles_octant, 0.0)
  , d_xi(d_number_angles_octant, 0.0)
  , d_name(name)
  , d_octant_sign(8, vec_dbl(3, 0.0))
{
  Insist(order > 0,
    "The quadrature order MUST be positive.");
  Insist((dim > 0) and (dim < 4),
    "The quadrature dimension MUST be 1, 2, or 3.");
  // Define the signs for all eight octants.
  // first
  d_octant_sign[0][0] =  1.0;
  d_octant_sign[0][1] =  1.0;
  d_octant_sign[0][2] =  1.0;
  // second
  d_octant_sign[1][0] = -1.0;
  d_octant_sign[1][1] =  1.0;
  d_octant_sign[1][2] =  1.0;
  // third
  d_octant_sign[2][0] = -1.0;
  d_octant_sign[2][1] = -1.0;
  d_octant_sign[2][2] =  1.0;
  // fourth
  d_octant_sign[3][0] =  1.0;
  d_octant_sign[3][1] = -1.0;
  d_octant_sign[3][2] =  1.0;
  // fifth
  d_octant_sign[4][0] =  1.0;
  d_octant_sign[4][1] =  1.0;
  d_octant_sign[4][2] = -1.0;
  // sixth
  d_octant_sign[5][0] = -1.0;
  d_octant_sign[5][1] =  1.0;
  d_octant_sign[5][2] = -1.0;
  // seventh
  d_octant_sign[6][0] = -1.0;
  d_octant_sign[6][1] = -1.0;
  d_octant_sign[6][2] = -1.0;
  // eighth
  d_octant_sign[7][0] =  1.0;
  d_octant_sign[7][1] = -1.0;
  d_octant_sign[7][2] = -1.0;
}

// Pure virtual still needs definition.
Quadrature::~Quadrature()
{
  /* ... */
}

// Display
void Quadrature::display() const
{

    using std::cout;
    using std::endl;
    using std::printf;

    cout << endl;
    cout << d_name << " abscissa and weights: " << endl << endl;

    double weight_sum = 0;

    if (d_dimension == 1)
    {
      cout << "   m            mu                  wt       " << endl;
      cout << "  ---   ------------------  -----------------" << endl;
      for ( int ix = 0; ix < d_number_angles_octant; ++ix )
      {
          printf ("%4i    %16.13f   %16.13f   \n", ix, d_mu[ix], d_weight[ix] );
          weight_sum += d_weight[ix];
      }
    }
    else
    {
      cout << "   m            mu                 eta                xi                 wt       " << endl;
      cout << "  ---   -----------------  -----------------  -----------------  -----------------" << endl;
      for ( int ix = 0; ix < d_number_angles_octant; ++ix )
      {
          printf ("%4i    %16.13f   %16.13f   %16.13f   %16.13f   \n",
                  ix, d_mu[ix], d_eta[ix], d_xi[ix], d_weight[ix] );
          weight_sum += d_weight[ix];
      }
    }

    cout << endl << "  The sum of the weights is " << weight_sum << endl;
    cout << endl;

}

}

//---------------------------------------------------------------------------//
//              end of Quadrature.cc
//---------------------------------------------------------------------------//
