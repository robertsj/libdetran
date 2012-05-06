//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2D.t.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper 2D specialization.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER2D_T_HH_
#define SWEEPER2D_T_HH_

// Detran
#include "Equation_DD_2D.hh"
#include "Equation_SC_2D.hh"
#include "detran_config.h"

// System
#include <iostream>

namespace detran
{
using std::cout;
using std::endl;

template<class EQ>
inline void Sweeper<EQ, _2D>::setup()
{

  // Check to see if input has equation type. If not, default is dd.
//  std::string equation = "dd";
//  if (d_input->check("equation"))
//  {
//    equation = d_input->get<std::string>("equation");
//  }
//  // Set the equation.
//  if (equation == "dd")
//  {
//    d_equation(EQ(d_mesh, d_material, d_quadrature, d_update_psi));
//  }
//  else if (equation == "sc")
//  {
//    d_equation(Equation_SC_2D(d_mesh, d_material, d_quadrature, d_update_psi));
//  }
//  else
//  {
//    THROW("Unsupported equation type");
//  }

  // Set up face/octant map.
  d_face_index.resize(4, vec2_int(2, vec_int(2, 0)));

  int inc[4][2] = {0, 2, 1, 2, 1, 3, 0, 3};
  int out[4][2] = {1, 3, 0, 3, 0, 2, 1, 2};
  for (int i = 0; i < 4; i++)
  {
    //      octant   surface type   inout
    d_face_index[i][Mesh::VERT][Boundary_T::IN]  = inc[i][0];
    d_face_index[i][Mesh::HORZ][Boundary_T::IN]  = inc[i][1];
    d_face_index[i][Mesh::VERT][Boundary_T::OUT] = out[i][0];
    d_face_index[i][Mesh::HORZ][Boundary_T::OUT] = out[i][1];
  }
}

#define STATIC_EQ


// Instantiate
template class Sweeper<Equation_DD_2D, _2D>;
template class Sweeper<Equation_SC_2D, _2D>;

} // end namespace detran

#endif /* SWEEPER2D_T_HH_ */
