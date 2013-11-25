//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CoarseMesh.cc
 *  @brief CoarseMesh member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "transport/CoarseMesh.hh"
#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "utilities/SoftEquivalence.hh"

namespace detran
{

//----------------------------------------------------------------------------//
CoarseMesh::CoarseMesh(SP_mesh fine_mesh, const size_t level)
  : d_fine_mesh(fine_mesh)
  , d_level(level)
  , d_fine_to_coarse(3)
  , d_coarse_edge_flag(3)
{
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

  d_fine_to_coarse[0].resize(1, 0);
  d_fine_to_coarse[1].resize(1, 0);
  d_fine_to_coarse[2].resize(1, 0);

  for (int d = 0; d < dim; d++)
  {

    // Try dividing the mesh by the level, recording any remainder.
    number_fine[d]   = d_fine_mesh->number_cells(d); // 5, level=3
    remainder[d]     = number_fine[d] % d_level;     // 5 % 3 = 2
    number_coarse[d] = (number_fine[d] - remainder[d]) / d_level;

    // Set fine meshes per coarse to be level.  If the meshing yields a single
    // coarse mesh, then the remainder is added to that.
    number_fine_coarse[d].resize(number_coarse[d], d_level);
    if (number_coarse[d] == 1)
    {
      number_fine_coarse[d][0] += remainder[d];
    }
    else
    {
      // Distribute any remaining fine meshes.
      for (int i = 0; i < remainder[d]; i++)
      {
        number_fine_coarse[d][i] += 1;
      }
    }

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
    d_coarse_mesh =
      new detran_geometry::Mesh1D(coarse_edges[0], coarse_material_map);
  }
  else if (dim == 2)
  {
    d_coarse_mesh =
      new detran_geometry::Mesh2D(coarse_edges[0], coarse_edges[1],
                                  coarse_material_map);
  }
  else
  {
    d_coarse_mesh =
      new detran_geometry::Mesh3D(coarse_edges[0], coarse_edges[1],
                                  coarse_edges[2], coarse_material_map);
  }

  // Create fine mesh to coarse mesh map
  vec_int f2c_mesh_map(d_fine_mesh->number_cells(), 0);
  for (size_t k = 0; k < d_fine_mesh->number_cells_z(); ++k)
  {
    size_t kk = this->fine_to_coarse(k, 2);
    for (size_t j = 0; j < d_fine_mesh->number_cells_y(); ++j)
    {
      size_t jj = this->fine_to_coarse(j, 1);
      for (size_t i = 0; i < d_fine_mesh->number_cells_x(); ++i)
      {
        size_t ii = this->fine_to_coarse(i, 0);
        f2c_mesh_map[d_fine_mesh->index(i, j, k)] =
          d_coarse_mesh->index(ii, jj, kk);
      }
    }
  }
  // \todo it might be good to check if the key is there and
  // use a numbered version
  d_fine_mesh->add_mesh_map("COARSEMESH", f2c_mesh_map);

  using detran_utilities::soft_equiv;
  Ensure(soft_equiv(d_fine_mesh->total_width_x(),
                    d_coarse_mesh->total_width_x()));
  Ensure(soft_equiv(d_fine_mesh->total_width_y(),
                    d_coarse_mesh->total_width_y()));
  Ensure(soft_equiv(d_fine_mesh->total_width_z(),
                    d_coarse_mesh->total_width_z()));
}

TRANSPORT_TEMPLATE_EXPORT(detran_utilities::SP<CoarseMesh>)

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file CoarseMesh.cc
//----------------------------------------------------------------------------//
