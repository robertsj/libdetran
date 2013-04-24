//----------------------------------*-C++-*----------------------------------//
/**
 * @file   Sweeper.t.hh
 * @author Jeremy Roberts
 * @date   Apr 22, 2012
 * @brief  Sweeper template member definitions
 */
//---------------------------------------------------------------------------//

#ifndef detran_SWEEPER_T_HH_
#define detran_SWEEPER_T_HH_

#include "transport/Sweeper.hh"

namespace detran
{

//---------------------------------------------------------------------------//
// Default setup
template <class D>
inline void Sweeper<D>::setup()
{
  THROW("NOT IMPLEMENTED");
}

//---------------------------------------------------------------------------//
// 3D Setup
template <>
inline void Sweeper<_3D>::setup()
{
  // Set up face/octant map.
  d_face_index.resize(8, vec2_size_t(3, vec_size_t(2, 0)));
  int inc[8][3] =
  {{0, 2, 4}, {1, 2, 4}, {1, 3, 4}, {0, 3, 4},
   {0, 2, 5}, {1, 2, 5}, {1, 3, 5}, {0, 3, 5}};
  int out[8][3] =
  {{1, 3, 5}, {0, 3, 5}, {0, 2, 5}, {1, 2, 5},
   {1, 3, 4}, {0, 3, 4}, {0, 2, 4}, {1, 2, 4}};
  for (int i = 0; i < 8; i++)
  {
    // octant     surface type    inout
    d_face_index[i][Mesh::YZ][Boundary_T::IN] = inc[i][0];
    d_face_index[i][Mesh::XZ][Boundary_T::IN] = inc[i][1];
    d_face_index[i][Mesh::XY][Boundary_T::IN] = inc[i][2];
    d_face_index[i][Mesh::YZ][Boundary_T::OUT] = out[i][0];
    d_face_index[i][Mesh::XZ][Boundary_T::OUT] = out[i][1];
    d_face_index[i][Mesh::XY][Boundary_T::OUT] = out[i][2];
  }
}

//---------------------------------------------------------------------------//
// 2D Setup
template <>
inline void Sweeper<_2D>::setup()
{
  // Set up face/octant map.
  d_face_index.resize(4, vec2_size_t(2, vec_size_t(2, 0)));
  // Example:  octant 1 is incident on faces 0 and 2
  //           and leaves faces 1 and 3.
  int inc[4][2] = {{0, 2}, {1, 2}, {1, 3}, {0, 3}};
  int out[4][2] = {{1, 3}, {0, 3}, {0, 2}, {1, 2}};
  for (int i = 0; i < 4; i++)
  {
    //      octant   surface type   inout
    d_face_index[i][Mesh::VERT][Boundary_T::IN]  = inc[i][0];
    d_face_index[i][Mesh::HORZ][Boundary_T::IN]  = inc[i][1];
    d_face_index[i][Mesh::VERT][Boundary_T::OUT] = out[i][0];
    d_face_index[i][Mesh::HORZ][Boundary_T::OUT] = out[i][1];
  }
}

//---------------------------------------------------------------------------//
// 1D Setup
template <>
inline void Sweeper<_1D>::setup()
{
  // Set up face/octant map.
  d_face_index.resize(2, vec2_size_t(1, vec_size_t(2, 0)));
  // octant / surface type / inout
  d_face_index[0][Mesh::VERT][Boundary_T::IN]  = Mesh::WEST;
  d_face_index[1][Mesh::VERT][Boundary_T::IN]  = Mesh::EAST;
  d_face_index[0][Mesh::VERT][Boundary_T::OUT] = Mesh::EAST;
  d_face_index[1][Mesh::VERT][Boundary_T::OUT] = Mesh::WEST;
}

} // end namespace detran

#endif /* detran_SWEEPER_T_HH_ */
