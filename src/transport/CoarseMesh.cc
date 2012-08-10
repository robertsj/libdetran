//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CoarseMesh.cc
 * \brief  CoarseMesh 
 * \author Jeremy Roberts
 * \date   Aug 8, 2012
 */
//---------------------------------------------------------------------------//

// Detran
#include "CoarseMesh.hh"
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

namespace detran
{

CoarseMesh::CoarseMesh(SP_mesh fine_mesh, u_int level)
  : d_fine_mesh(fine_mesh)
  , d_level(level)
  , d_fine_to_coarse(3)
  , d_coarse_edge_flag(3)
{
  // Preconditions
  Require(d_fine_mesh);

  using std::cout;
  using std::endl;

  // Problem dimension
  int dim = d_fine_mesh->dimension();

  // Number of fine meshes per coarse mesh
  vec_int number_fine(dim);
  vec_int remainder(dim);
  vec_int number_coarse(dim);
  vec2_int number_fine_coarse(dim);
  vec2_dbl coarse_edges(dim);

  d_coarse_edge_flag[0].resize(d_fine_mesh->number_cells_x() + 1, -1);
  d_coarse_edge_flag[1].resize(d_fine_mesh->number_cells_y() + 1, -1);
  d_coarse_edge_flag[2].resize(d_fine_mesh->number_cells_z() + 1, -1);

  for (int d = 0; d < dim; d++)
  {

    // Try dividing the mesh by the level, recording any remainder.
    number_fine[d]   = d_fine_mesh->number_cells(d);
    remainder[d]     = number_fine[d] % level;
    number_coarse[d] = (number_fine[d] - remainder[d]) / level;

    // Set fine meshes per coarse to be level.
    number_fine_coarse[d].resize(number_coarse[d], level);

    // Distribute any remaining fine meshes.
    for (int i = 0; i < remainder[d]; i++)
      number_fine_coarse[d][i] += 1;

    // Resize the edge and fine to coarse vectors.
    coarse_edges[d].resize(number_coarse[d] + 1, 0.0);
    d_fine_to_coarse[d].resize(number_fine[d], 0);

    // Set the first fine mesh edge to coincide with the first coarse mesh edge.
    int edge_flag = 0;
    d_coarse_edge_flag[d][0] = edge_flag++;

    // Fill the mesh edges and the fine to coarse vectors.
    int f_start = 0;
    int f_stop = 0;
    for (int c = 0; c < number_coarse[d]; c++)
    {
      double edge = coarse_edges[d][c];
      f_stop += number_fine_coarse[d][c];
      for (int f = f_start; f < f_stop; f++)
      {
        edge += d_fine_mesh->width(d, f);
        d_fine_to_coarse[d][f] = c;
      }
      coarse_edges[d][c + 1] = edge;
      f_start += number_fine_coarse[d][c];
      // Set the fine mesh edge to match a coarse mesh edge.
      Assert(f_stop < d_coarse_edge_flag[d].size());
      d_coarse_edge_flag[d][f_stop] = edge_flag++;
    }
    Assert(edge_flag == coarse_edges[d].size());

  } // end dimension loop

  // Create the mesh.
  int number_coarse_mesh = 1;
  for (int d = 0; d < dim; d++)
    number_coarse_mesh *= number_coarse[d];

  // Each coarse mesh material is unique.
  vec_int coarse_material_map(number_coarse_mesh, 0);
  for (int c = 0; c < number_coarse_mesh; c++)
    coarse_material_map[c] = c;

  // Create the coarse mesh.
  if (dim == 1)
  {
    d_coarse_mesh = new Mesh1D(coarse_edges[0],
                               coarse_material_map);
  }
  else if (dim == 2)
  {
    d_coarse_mesh = new Mesh2D(coarse_edges[0],
                               coarse_edges[1],
                               coarse_material_map);
  }
  else
  {
    d_coarse_mesh = new Mesh3D(coarse_edges[0],
                               coarse_edges[1],
                               coarse_edges[2],
                               coarse_material_map);
  }

}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file CoarseMesh.cc
//---------------------------------------------------------------------------//
