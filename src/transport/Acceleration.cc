//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Acceleration.cc
 * \author Jeremy Roberts
 * \date   May 17, 2012
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
Acceleration::Acceleration(SP_mesh mesh,
                           SP_material material,
                           SP_quadrature quadrature)
  : b_mesh(mesh)
  , b_material(material)
  , b_quadrature(quadrature)
{
  Require(b_mesh);
  Require(b_material);

}

// Create coarse mesh.  Note, I'm trying to avoid any special
// treatment for dimension.
// \todo I already see room for refactoring.  Do xyz as loop.
void Acceleration::initialize(int level)
{
  using std::cout;
  using std::endl;

  Require(level > 0);

  int number_fine_x = b_mesh->number_cells_x();
  int number_fine_y = b_mesh->number_cells_y();
  int number_fine_z = b_mesh->number_cells_z();

  // Compute first try of the number of coarse meshes, keeping
  // the remainder for round-robin allocation
  int remainder_x = number_fine_x % level;
  int number_coarse_x = (number_fine_x - remainder_x) / level;
  cout << " remainder x = " << remainder_x
       << " number_coarse_x " << number_coarse_x
       << " number_fine_x " << number_fine_x
       << " level = " << level << endl;
//  int remainder_y = number_fine_y % level;
//  int number_coarse_y = (number_fine_y - remainder_y) / level;
//  int remainder_z = number_fine_z % level;
//  int number_coarse_z = (number_fine_z - remainder_z) / level;
  int number_coarse_y = 1;
  int number_coarse_z = 1;
  // Initial number of fine per coarse mesh in each direction
  vec_int number_fine_coarse_x(number_coarse_x, level);
//  vec_int number_fine_coarse_y(number_coarse_y, level);
//  vec_int number_fine_coarse_z(number_coarse_z, level);

  // Add the remainders.
  for (int i = 0; i < remainder_x; i++)
  {
    number_fine_coarse_x[i] += 1;
  }
//  for (int i = 0; i < remainder_y; i++)
//  {
//    number_fine_coarse_y[i] += 1;
//  }
//  for (int i = 0; i < remainder_z; i++)
//  {
//    number_fine_coarse_z[i] += 1;
//  }

  // Create the edges
  vec_dbl coarse_edges_x(number_coarse_x + 1, 0.0);
//  vec_dbl coarse_edges_y(number_coarse_y + 1, 0.0);
//  vec_dbl coarse_edges_z(number_coarse_z + 1, 0.0);
  // x
  int f_holder = 0;
  int f_stop = 0;
  for (int c = 0; c < number_coarse_x; c++)
  {
    double edge = coarse_edges_x[c];
    f_stop += number_fine_coarse_x[c];
    cout << "f_holder = " << f_holder << endl;
    for (int f = f_holder; f < f_stop; f++)
    {
      cout << "*f_holder = " << f_holder << endl;
      edge += b_mesh->dx(f);
      cout << " edge = " << edge << endl;
    }
    coarse_edges_x[c + 1] = edge;
    f_holder += number_fine_coarse_x[c];
  }
//  // y
//  f_holder = 0;
//  for (int c = 0; c < number_coarse_y; c++)
//  {
//    double edge = coarse_edges_y[c];
//    for (int f = f_holder; f < number_fine_coarse_y[c]; f++)
//    {
//      edge += b_mesh->dy(f);
//    }
//    coarse_edges_y[c + 1] = edge;
//    f_holder += number_fine_coarse_y[c];
//  }
//  // z
//  f_holder = 0;
//  for (int c = 0; c < number_coarse_z; c++)
//  {
//    double edge = coarse_edges_z[c];
//    for (int f = f_holder; f < number_fine_coarse_z[c]; f++)
//    {
//      edge += b_mesh->dz(f);
//    }
//    coarse_edges_z[c + 1] = edge;
//    f_holder += number_fine_coarse_z[c];
//  }

  // Coarse mesh materials.  Assume for now that the number of
  // groups does not change.
  int number_coarse_mesh = number_coarse_x * number_coarse_y * number_coarse_z;
  b_coarse_material =
      new Material(number_coarse_mesh, b_material->number_groups());

  // Each coarse mesh material is unique.
  vec_int coarse_material_map(number_coarse_mesh, 0);
  for (int c = 0; c < number_coarse_mesh; c++)
  {
    coarse_material_map[c] = c;
  }

  if (b_mesh->dimension() == 1)
  {
     b_coarse_mesh = new Mesh1D(coarse_edges_x,
                                coarse_material_map);
  }
//  else if (b_mesh->dimension() == 2)
//  {
//    b_coarse_mesh = new Mesh2D(coarse_edges_x,
//                               coarse_edges_y,
//                               coarse_material_map);
//  }
//  else
//  {
//    b_coarse_mesh = new Mesh3D(coarse_edges_x,
//                               coarse_edges_y,
//                               coarse_edges_z,
//                               coarse_material_map);
//  }


  // Make the appropriate mesh.


}


void Acceleration::homogenize(SP_state state)
{
  Require(state);

  // Loop through all coarse mesh cells
  //   Loop through all fine mesh cells in the coarse mesh cell
  //     Loop through all groups
  //       Compute the average group constants

}

} // end namespace detran


