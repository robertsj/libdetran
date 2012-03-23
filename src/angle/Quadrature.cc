//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Quadrature.cc
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//

// Angle headers
#include "Quadrature.hh"

// System headers
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
  Require(dim > 0);
  Require(dim < 4);
  Require(order > 0);
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

}

//---------------------------------------------------------------------------//
//              end of Quadrature.cc
//---------------------------------------------------------------------------//
