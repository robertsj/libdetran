//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper1D.t.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper 1D specialization.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//
#ifndef SWEEPER1D_T_HH_
#define SWEEPER1D_T_HH_

namespace detran
{

// Constructor setup.
template <>
inline void Sweeper<_1D>::setup(SP_material material)
{
  // Set up face/octant map.
  d_face_index.resize(2, vec2_int(1, vec_int(2, 0)));
  // octant / surface type / inout
  d_face_index[0][Mesh::VERT][Boundary_T::IN] = Mesh::LEFT;
  d_face_index[1][Mesh::VERT][Boundary_T::IN] = Mesh::RIGHT;
  d_face_index[0][Mesh::VERT][Boundary_T::OUT] = Mesh::RIGHT;
  d_face_index[1][Mesh::VERT][Boundary_T::OUT] = Mesh::LEFT;
}

// Instantiate
template class Sweeper<_1D> ;

} // end namespace detran

#endif /* SWEEPER1D_T_HH_ */
