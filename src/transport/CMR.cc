/*
 * CMR.cc
 *
 *  Created on: May 17, 2012
 *      Author: robertsj
 */

// Detran
#include "CMR.hh"

// System
#include <algorithm>

namespace detran
{

template <class D>
CMR<D>::CMR(SP_input input,
            SP_material material,
            SP_coarsemesh coarsemesh,
            SP_currenttally currenttally)
  : Base(input, material, coarsemesh, currenttally)
{

}

template <class D>
void CMR<D>::update(State::moments_type &phi, SP_sweepsource source)
{
  Require(source);

  // Integrate the coarse meshes
  integrate(phi, source);

  return;
}

//template <>
//void CMR<_1D>::update(State::moments_type &phi, SP_sweepsource source)
//{
//  return;
//}

//
//    Require(source);
//
//    // Save the current fine mesh flux.
//    State::moments_type phi_old(phi);
//
//    // Integrate the coarse meshes
//    integrate(phi, source);
//
//    // Compute the linear system (template specialization)
//
//    // Solve the linear system for the rebalance factors
//
//    // Update the flux
//
//}

//----------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//----------------------------------------------------------------------------//

template <class D>
void CMR<D>::integrate(State::moments_type &phi, SP_sweepsource source)
{
  int number_coarse = b_coarse_mesh->number_cells();

   //vec_dbl phi_bar(number_coarse, 0.0);
   //vec_dbl sigma_t(number_coarse, 0.0);

   for (int cell_coarse = 0; cell_coarse < number_coarse; cell_coarse++)
   {

     // Temporary homogenized constants.
     double phi_bar = 0.0;
     double volume_coarse = 0.0;

     for (int k = 0; k < b_mesh->number_cells_z(); k++)
     {
       for (int j = 0; j < b_mesh->number_cells_y(); j++)
       {
         for (int i = 0; i < b_mesh->number_cells_x(); i++)
         {
           int cell_fine = b_mesh->index(i, j, k);
           double volume_fine = b_mesh->dx(i) * b_mesh->dy(j) * b_mesh->dz(k);
           volume_coarse += volume_fine;
           phi_bar += volume_fine * phi[cell_fine];
         }
       }
     }
     phi_bar /= volume_coarse;

   } // end coarse cell loop
}

template <>
void CMR<_1D>::integrate(State::moments_type &phi, SP_sweepsource sweepsource)
{
  int number_coarse = b_coarse_mesh->number_cells();

  vec_int mat_map = b_mesh->mesh_map("MATERIAL");

  for (int cell_coarse = 0; cell_coarse < number_coarse; cell_coarse++)
  {

    double vol_c = 0.0;

    for (int i = 0; i < b_mesh->number_cells_x(); i++)
    {

      // Fine mesh index.
      int cell = i;

      // Fine mesh material index.
      int mat = mat_map[cell];

      // Fine mesh volume and contribution to coarse mesh volume.
      double vol_f = b_mesh->dx(i);
      vol_c += vol_f;

      // Add the fine mesh removal reaction rate, equal to the
      // total reaction rate minus the within-group scattering rate.
      d_R[cell_coarse] += vol_f * phi[cell] *
          (b_material->sigma_t(mat, b_g) - b_material->sigma_s(mat, b_g, b_g));


      // Remove the within-group scattering source from the right hand side.
      d_Q[cell_coarse] -= phi[cell] * b_material->sigma_s(mat, b_g, b_g);

      // Add up the discrete angular sources.
      for (int o = 0; o < 2; o++)
      {
        for (int a = 0; a < b_quadrature->number_angles(); a++)
        {
          // Get the angle source
          SweepSource<_1D>::sweep_source_type source =
              sweepsource->source(b_g, o, a);
          d_Q[cell_coarse] += b_quadrature->weight(a) * source[cell];
        }
      }

      d_Q[cell_coarse] *= vol_f;

    }


  } // end coarse cell loop

}

template <class D>
void CMR<D>::reset()
{
  std::fill(d_J_pos.begin(), d_J_pos.end(), 0.0);
  std::fill(d_J_neg.begin(), d_J_neg.end(), 0.0);
}

} // end namespace detran


