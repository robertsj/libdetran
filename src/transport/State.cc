//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   State.cc
 *  @author Jeremy Roberts
 *  @date   Mar 25, 2012
 *  @brief  State member definitions.
 */
//---------------------------------------------------------------------------//

#include "State.hh"
#include <iostream>
#include <cstdio>

namespace detran
{

//---------------------------------------------------------------------------//
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
  Require(d_number_groups > 0);

  // Allocate angular flux vectors if needed.
  int store_psi = 0;
  if (input->check("store_angular_flux"))
  {
    store_psi = input->get<int>("store_angular_flux");
  }

  if (store_psi > 0)
  {
    Insist(d_quadrature, "Angular flux requested but no quadrature given.");

    d_store_angular_flux = true;

    for (int g = 0; g < d_number_groups; g++)
    {
      d_angular_flux[g].resize(d_quadrature->number_angles(),
                               vec_dbl(mesh->number_cells(), 0.0));
    }

  }

}

//---------------------------------------------------------------------------//
void State::clear()
{
  for (int g = 0; g < d_number_groups; ++g)
  {
    for (int i = 0; i < d_mesh->number_cells(); ++i)
    {
      d_moments[g][i] = 0.0;
      if (d_store_angular_flux)
      {
        for (int a = 0; a < d_quadrature->number_angles(); ++a)
        {
          d_angular_flux[g][a][i] = 0.0;
        }
      }
    }
  }

}

//---------------------------------------------------------------------------//
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

  printf("\n");
  printf("\n");
  for (int g = 0; g < d_moments.size() + 1; g++)
    printf("--------------");
  printf("\n");
  printf("Scalar Flux Moments\n");
  printf("cell \\ g");
  for (int g = 0; g < d_number_groups; g++)
    printf(" %12i ", g);
  printf("\n");
  for (int g = 0; g < d_moments.size() + 1; g++)
    printf("--------------");
  printf("\n");
  for (int i = 0; i < d_mesh->number_cells(); i++)
  {
    printf("%10i", i);
    for (int g = 0; g < d_number_groups; g++)
    {
      printf(" %12.5e ", d_moments[g][i]);
    }
    printf("\n");
  }

  if (d_store_angular_flux)
  {
    printf("\n");
    for (int a = 0; a < d_angular_flux[0].size() + 1; a++)
      printf("--------------", a);
    printf("\n");
    printf("Discrete Angular Flux\n");
    for (int a = 0; a < d_angular_flux[0].size() + 1; a++)
      printf("--------------", a);
    printf("\n");

    for (int g = 0; g < d_number_groups; g++)
    {
      printf("group %4i \n", g);
      printf("cell \\ a");
      for (int a = 0; a < d_angular_flux[g].size(); a++)
        printf(" %12i ", a);
      printf("\n");
      for (int a = 0; a < d_angular_flux[0].size() + 1; a++)
        printf("--------------", a);
      printf("\n");
      for (int i = 0; i < d_mesh->number_cells(); i++)
      {
        printf("%10i", i);
        for (int a = 0; a < d_angular_flux[g].size(); a++)
        {
          printf(" %12.5e ", d_angular_flux[g][a][i]);
        }
        printf("\n");
      }
      printf("\n");
    }
  }

}

}

//---------------------------------------------------------------------------//
//              end of State.cc
//---------------------------------------------------------------------------//
