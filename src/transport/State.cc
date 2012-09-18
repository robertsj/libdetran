//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   State.cc
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  State member definitions.
 */
//---------------------------------------------------------------------------//

#include "State.hh"
#include <iostream>
#include <cstdio>

namespace detran
{

// Constructor.
State::State(SP_input        input,
             SP_mesh         mesh,
             SP_quadrature   quadrature)
  : d_input(input)
  , d_mesh(mesh)
  , d_quadrature(quadrature)
  , d_number_groups(input->get<int>("number_groups"))
  , d_store_angular_flux(false)
  , d_moments(d_number_groups, vec_dbl(mesh->number_cells(), 0.0))
  , d_angular_flux(d_number_groups)
  , d_eigenvalue(0.0)
{
  // Preconditions
  Require(input);
  Require(mesh);
  Require(quadrature);
  Require(d_number_groups > 0);

  // Allocate angular flux vectors if needed.
  int store_psi = 0;
  if (input->check("store_angular_flux"))
  {
    store_psi = input->get<int>("store_angular_flux");
  }

  if (store_psi > 0)
  {
    d_store_angular_flux = true;

    for (int g = 0; g < d_number_groups; g++)
    {
      d_angular_flux[g].resize(quadrature->number_angles(),
                               vec_dbl(mesh->number_cells(), 0.0));
    }

  }

}

// Constructor.
State::State(SP_input        input,
             SP_mesh         mesh)
  : d_input(input)
  , d_mesh(mesh)
  , d_number_groups(input->get<int>("number_groups"))
  , d_store_angular_flux(false)
  , d_moments(d_number_groups, vec_dbl(mesh->number_cells(), 0.0))
  , d_eigenvalue(0.0)
{
  // Preconditions
  Require(input);
  Require(mesh);
  Require(d_number_groups > 0);
}

void State::display() const
{
  using std::printf;
  printf("\n");
  printf("------------\n");
  printf("Detran State\n");
  printf("------------\n");
  printf(" number groups: %5i \n", d_number_groups);
  printf("  number cells: %5i \n", d_mesh->number_cells());
  printf("     store psi: %5i \n", d_store_angular_flux);
  // print phi as function of cell and group
  printf("\n------------------------------------------------------------\n");
  printf(" cell \\ g    ");
  for (int g = 0; g < d_number_groups; g++)
    printf(" %4i ", g);
  printf("\n------------------------------------------------------------\n");
  for (int i = 0; i < d_mesh->number_cells(); i++)
  {
    printf(" %5i ", i);
    for (int g = 0; g < d_number_groups; g++)
    {
      printf(" %12.5e ", d_moments[g][i]);
    }
    printf("\n");
  }
}

}

//---------------------------------------------------------------------------//
//              end of State.cc
//---------------------------------------------------------------------------//
