//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Reflective.i.hh
 * \author robertsj
 * \date   Apr 10, 2012
 * \brief  Reflective inline member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef REFLECTIVE_I_HH_
#define REFLECTIVE_I_HH_

// System
#include <iostream>

namespace detran
{

template <class D>
void Reflective<D>::update(int g)
{
  for (int o = 0; o < d_quadrature->number_octants()/2; o++)
  {
    for (int a = 0; a < d_quadrature->number_angles_octant(); a++)
    {
      // Reroute reflecting fluxes.
      d_boundary(d_side, d_octants[o][Boundary_T::IN], a, g) =
        d_boundary(d_side, d_octants[o][Boundary_T::OUT], a, g);
    }
  }
}

template <class D>
void Reflective<D>::update(int g, int o, int a)
{
  int o_out = -1;
  for (int i = 0; i < d_octants.size(); i++)
  {
    if (d_octants[i][Boundary_T::IN] == o)
    {
      o_out = d_octants[i][Boundary_T::OUT];
    }
  }
  // Only reroute fluxes if I am an outgoing side
  if (o_out >= 0)
  {
    // Reroute reflecting fluxes.
    d_boundary(d_side, o, a, g) =
      d_boundary(d_side, o_out, a, g);
  }
}

//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//

// Note, all octants are arrange left-to-right-up, left-to-right-down
// with respect to the incident normal.  Left-to-right is defined by
// the secondary axis (x->y->z, y->z->x, z->x->y) and the normal
// right hand rule.  For example, if I'm
// incident on the left face, I am looking along the +x, my left is
// oriented along +y, and above me is +z.  If I'm incident on the south
// face, I'm looking along +z, my left is -x
template <class D>
void Reflective<D>::setup_octant()
{
  if (d_side == Mesh::LEFT)
  {
    // incident octants: 0, 3, 4, 7
    d_octants[0][0] = 0;
    d_octants[1][0] = 3;
    d_octants[2][0] = 4;
    d_octants[3][0] = 7;
    // outgoing octants: 1, 2, 5, 6
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
    d_octants[2][1] = 5;
    d_octants[3][1] = 6;
  }
  else if (d_side == Mesh::RIGHT)
  {
    // incident octants: 2, 1, 6, 5
    d_octants[0][0] = 2;
    d_octants[1][0] = 1;
    d_octants[2][0] = 6;
    d_octants[3][0] = 5;
    // outgoing octants: 0, 3, 4, 7
    d_octants[0][1] = 0;
    d_octants[1][1] = 3;
    d_octants[2][1] = 4;
    d_octants[3][1] = 7;
  }
  else if (d_side == Mesh::BOTTOM)
  {
    // incident octants: 0, 4, 1, 5
    d_octants[0][0] = 0;
    d_octants[1][0] = 4;
    d_octants[2][0] = 1;
    d_octants[3][0] = 5;
    // outgoing octants: 3, 7, 2, 6
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
    d_octants[2][1] = 5;
    d_octants[3][1] = 6;
  }
  else if (d_side == Mesh::TOP)
  {
    // incident octants: 7, 3, 6, 2
    d_octants[0][0] = 0;
    d_octants[1][0] = 3;
    d_octants[2][0] = 4;
    d_octants[3][0] = 7;
    // outgoing octants: 4, 0, 5, 1
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
    d_octants[2][1] = 5;
    d_octants[3][1] = 6;
  }
  else if (d_side == Mesh::SOUTH)
  {
    // incident octants: 1, 0, 2, 3
    d_octants[0][0] = 0;
    d_octants[1][0] = 3;
    d_octants[2][0] = 4;
    d_octants[3][0] = 7;
    // outgoing octants: 5, 4, 6, 7
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
    d_octants[2][1] = 5;
    d_octants[3][1] = 6;
  }
  else
  {
    // incident octants: 4, 5, 7, 6
    d_octants[0][0] = 0;
    d_octants[1][0] = 3;
    d_octants[2][0] = 4;
    d_octants[3][0] = 7;
    // outgoing octants: 3, 2, 0, 1
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
    d_octants[2][1] = 5;
    d_octants[3][1] = 6;
  }
}

template <>
void Reflective<_2D>::setup_octant()
{
  Require(d_side == Mesh::LEFT   || d_side == Mesh::RIGHT ||
          d_side == Mesh::BOTTOM || d_side == Mesh::TOP   );
  if (d_side == Mesh::LEFT)
  {
    d_octants[0][0] = 0; // [left/right oct][in/out]
    d_octants[1][0] = 3;
    d_octants[0][1] = 1;
    d_octants[1][1] = 2;
  }
  else if (d_side == Mesh::RIGHT)
  {
    d_octants[0][0] = 2;
    d_octants[1][0] = 1;
    d_octants[0][1] = 3;
    d_octants[1][1] = 0;
  }
  else if (d_side == Mesh::BOTTOM)
  {
    d_octants[0][0] = 1;
    d_octants[1][0] = 0;
    d_octants[0][1] = 2;
    d_octants[1][1] = 3;
  }
  else // TOP
  {
    d_octants[0][0] = 3;
    d_octants[1][0] = 2;
    d_octants[0][1] = 0;
    d_octants[1][1] = 1;
  }
}


template <>
void Reflective<_1D>::setup_octant()
{
  Require(d_side == Mesh::LEFT || d_side == Mesh::RIGHT);
  if (d_side == Mesh::LEFT)
  {
    // incident octant: 0
    d_octants[0][0] = 0;
    // outgoing octant: 1
    d_octants[0][1] = 1;
  }
  else
  {
    // incident octant: 1
    d_octants[0][0] = 1;
    // outgoing octant: 0
    d_octants[0][1] = 0;
  }
}

} // end namespace detran

#endif /* REFLECTIVE_I_HH_ */
