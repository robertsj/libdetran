/*
 * Sweeper.t.hh
 *
 *  Created on: Apr 22, 2012
 *      Author: robertsj
 */

#ifndef SWEEPER_T_HH_
#define SWEEPER_T_HH_

namespace detran
{

// Constructor
template <class D>
Sweeper<D>::Sweeper(SP_input input,
                    SP_mesh mesh,
                    SP_material material,
                    SP_quadrature quadrature,
                    SP_state state,
                    SP_sweepsource sweepsource)
  : d_input(input)
  , d_mesh(mesh)
  , d_material(material)
  , d_quadrature(quadrature)
  , d_state(state)
  , d_sweepsource(sweepsource)
  , d_update_psi(false)
  , d_adjoint(false)
  , d_number_sweeps(0)
  , d_update_boundary(false)
{
  Require(d_input);
  Require(d_mesh);
  Require(d_material);
  Require(d_quadrature);
  Require(d_state);
  Require(d_sweepsource);
  setup();
  // Check whether we keep psi.
  if (d_input->check("store_angular_flux"))
  {
    d_update_psi = d_input->get<int>("store_angular_flux");
  }
}

// 3D Setup
template <class D>
inline void Sweeper<D>::setup()
{

  // Set up face/octant map.
  d_face_index.resize(8, vec2_int(3, vec_int(2, 0)));

  int inc[8][3] =
  { 0, 2, 4, 1, 2, 4, 1, 3, 4, 0, 3, 4, 0, 2, 5, 1, 2, 5, 1, 3, 5, 0, 3, 5 };
  int out[8][3] =
  { 1, 3, 5, 0, 3, 5, 0, 2, 5, 1, 2, 5, 1, 3, 4, 0, 3, 4, 0, 2, 4, 1, 2, 4 };
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

// 2D Setup
template<>
inline void Sweeper<_2D>::setup()
{
  // Set up face/octant map.
  d_face_index.resize(4, vec2_int(2, vec_int(2, 0)));
  // Example:  octant 1 is incident on faces 0 and 2
  //           and leaves faces 1 and 3.
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

// 1D Setup
template <>
inline void Sweeper<_1D>::setup()
{
  // Set up face/octant map.
  d_face_index.resize(2, vec2_int(1, vec_int(2, 0)));
  // octant / surface type / inout
  d_face_index[0][Mesh::VERT][Boundary_T::IN] = Mesh::LEFT;
  d_face_index[1][Mesh::VERT][Boundary_T::IN] = Mesh::RIGHT;
  d_face_index[0][Mesh::VERT][Boundary_T::OUT] = Mesh::RIGHT;
  d_face_index[1][Mesh::VERT][Boundary_T::OUT] = Mesh::LEFT;
}

// Mesh sweeper indices
template <class D>
inline int Sweeper<D>::index(int o, int dim, int ijk)
{
  if (dim == 1)
  {
    if ((o == 0 or o == 3 or o == 4 or o == 7))
      return ijk;
    else
      return d_mesh->number_cells_x() - ijk - 1;
  }
  if (dim == 2)
  {
    if ((o == 0 or o == 1 or o == 4 or o == 5))
      return ijk;
    else
      return d_mesh->number_cells_y() - ijk - 1;
  }
  if (dim == 3)
  {
    if ((o == 0 or o == 1 or o == 2 or o == 3))
      return ijk;
    else
      return d_mesh->number_cells_z() - ijk - 1;
  }
}

} // end namespace detran

#endif /* SWEEPER_T_HH_ */
