//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Acceleration.cc
 * \author Jeremy Roberts
 * @date   May 17, 2012
 * \brief  Acceleration member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#include "Acceleration.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

#include <iostream>

namespace detran
{

// Constructor
template <class D>
Acceleration<D>::Acceleration(SP_mesh mesh,
                              SP_material material,
                              SP_quadrature quadrature)
  : b_mesh(mesh)
  , b_material(material)
  , b_quadrature(quadrature)
{
  Require(b_mesh);
  Require(b_material);
}

//----------------------------------------------------------------------------//
// Implementation
//----------------------------------------------------------------------------//

template <class D>
void Acceleration<D>::coarsen(int level)
{

  using std::cout;
  using std::endl;

  Require(level > 0);

  int dim = D::dimension;

  // Number of fine meshes per coarse mesh
  int number_fine[dim];
  int remainder[dim];
  int number_coarse[dim];
  vec2_int number_fine_coarse(dim);
  vec2_dbl coarse_edges(dim);

  b_fine_to_coarse.resize(dim);
  b_coarse_edge_flag.resize(dim);

  for (int d = 0; d < dim; d++)
  {
    // Try dividing the fine mesh by the desired level, recording
    // any remainder.
    number_fine[d]   = b_mesh->number_cells_x();
    remainder[d]     = number_fine[d] % level;
    number_coarse[d] = (number_fine[d] - remainder[d]) / level;

    // Set fine meshes per coarse to be level.
    number_fine_coarse[d].resize(number_coarse[d], level);

    // Distribute any remaining fine meshes.
    for (int i = 0; i < remainder[d]; i++)
    {
      number_fine_coarse[d][i] += 1;
    }

    // Resize the edge and fine to coarse vectors.
    coarse_edges[d].resize(number_coarse[d] + 1, 0.0);
    b_fine_to_coarse[d].resize(number_fine[d], 0);
    b_coarse_edge_flag[d].resize(number_fine[d], 0);

    // Fill the mesh edges and the fine to coarse vectors.
    int f_start = 0;
    int f_stop = 0;
    for (int c = 0; c < number_coarse[d]; c++)
    {
      double edge = coarse_edges[d][c];
      f_stop += number_fine_coarse[d][c];
      cout << "f_holder = " << f_start << endl;
      for (int f = f_start; f < f_stop; f++)
      {
        edge += b_mesh->width(d, f);
        cout << " edge = " << edge << endl;
        b_fine_to_coarse[d][f] = c;
      }
      coarse_edges[d][c + 1] = edge;
      f_start += number_fine_coarse[d][c];
    }
  }

  // Octant shift
  b_octant_shift.resize(3, vec_int(8, 0));
  // mu
  b_octant_shift[0][0] = 0;
  b_octant_shift[0][1] = 1;
  b_octant_shift[0][2] = 1;
  b_octant_shift[0][3] = 0;
  b_octant_shift[0][4] = 0;
  b_octant_shift[0][5] = 1;
  b_octant_shift[0][6] = 1;
  b_octant_shift[0][7] = 0;
  // eta
  b_octant_shift[1][0] = 0;
  b_octant_shift[1][1] = 0;
  b_octant_shift[1][2] = 1;
  b_octant_shift[1][3] = 1;
  b_octant_shift[1][4] = 0;
  b_octant_shift[1][5] = 0;
  b_octant_shift[1][6] = 1;
  b_octant_shift[1][7] = 1;
  // xi
  b_octant_shift[2][0] = 0;
  b_octant_shift[2][1] = 0;
  b_octant_shift[2][2] = 0;
  b_octant_shift[2][3] = 0;
  b_octant_shift[2][4] = 1;
  b_octant_shift[2][5] = 1;
  b_octant_shift[2][6] = 1;
  b_octant_shift[2][7] = 1;

  // Coarse mesh materials.  Assume for now that the number of
  // groups does not change.
  int number_coarse_mesh = 1;
  for (int d = 0; d < dim; d++)
  {
    number_coarse_mesh *= number_coarse[d];
  }
//  b_coarse_material =
//      new Material(number_coarse_mesh, b_material->number_groups());

  // Each coarse mesh material is unique.
  vec_int coarse_material_map(number_coarse_mesh, 0);
  for (int c = 0; c < number_coarse_mesh; c++)
  {
    coarse_material_map[c] = c;
  }

  // Create the coarse mesh.
  if (dim == 1)
  {
    b_coarse_mesh = new Mesh1D(coarse_edges[0],
                               coarse_material_map);
  }
  else if (dim == 2)
  {
    b_coarse_mesh = new Mesh2D(coarse_edges[0],
                               coarse_edges[1],
                               coarse_material_map);
  }
  else
  {
    b_coarse_mesh = new Mesh3D(coarse_edges[0],
                               coarse_edges[1],
                               coarse_edges[2],
                               coarse_material_map);
  }


}

//template <class D>
//void Acceleration<D>::homogenize(SP_state state, int group)
//{
//  Require(state);
//
//  // Get the group flux.
//  State::moments_type phi_g = state->phi(group);
//
//  int number_coarse = b_coarse_mesh->number_cells();
//
//  //vec_dbl phi_bar(number_coarse, 0.0);
//  //vec_dbl sigma_t(number_coarse, 0.0);
//
//  for (int cell_coarse = 0; cell_coarse < number_coarse; cell_coarse++)
//  {
//
//    // Temporary homogenized constants.
//    double phi_bar = 0.0;
//    double volume_coarse = 0.0;
//
//    for (int k = 0; k < b_mesh->number_cells_z(); k++)
//    {
//      for (int j = 0; j < b_mesh->number_cells_y(); j++)
//      {
//        for (int i = 0; i < b_mesh->number_cells_x(); i++)
//        {
//          int cell_fine = b_mesh->index(i, j, k);
//          double volume_fine = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);
//          volume_coarse += volume_fine;
//          phi_bar += volume_fine * phi_g[cell_fine];
//        }
//      }
//    }
//    phi_bar /= volume_coarse;
//
//  } // end coarse cell loop
//
//  // Loop through all coarse mesh cells
//  //   Loop through all fine mesh cells in the coarse mesh cell
//  //     Loop through all groups
//  //       Compute the average group constants
//
//}

} // end namespace detran


