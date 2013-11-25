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
  , d_number_groups(0)
  , d_number_moments(1)
  , d_eigenvalue(0.0)
  , d_store_angular_flux(false)
  , d_store_current(false)
  , d_adjoint(false)
{
  Require(d_input);
  Require(d_mesh);

  Insist(d_input->check("number_groups"),
         "Input must specify number of groups");
  d_number_groups = d_input->get<int>("number_groups");
  Require(d_number_groups > 0);

  // If a discrete problem, define the Legendre order
  if (d_quadrature)
  {
    size_t order = 0;
    if (input->check("moment_order"))
      order = input->get<int>("moment_order");
    d_momentindexer =
        detran_angle::MomentIndexer::Create(d_mesh->dimension(), order);
    d_number_moments = d_momentindexer->number_moments();
  }

  // Allocate the scalar flux moments vectors
  d_moments.resize(d_number_groups, 
                   vec_dbl(d_mesh->number_cells() * d_number_moments, 0.0));

  //  Allocate the current vectors
  int store_current = 0;
  if (input->check("store_current"))
    store_current = input->get<int>("store_current");
  if (store_current > 0)
  {
    d_store_current = true;
    d_current.resize(d_number_groups,
                     vec_dbl(d_mesh->number_cells(), 0.0));
  }

  // Allocate the angular flux vectors if needed.
  int store_psi = 0;
  if (input->check("store_angular_flux"))
    store_psi = input->get<int>("store_angular_flux");
  if (store_psi > 0)
  {
    Insist(d_quadrature, "Angular flux requested but no quadrature given.");
    d_store_angular_flux = true;
    d_angular_flux.resize(d_number_groups,
                          vec_moments_type(d_quadrature->number_angles(),
			               moments_type(d_mesh->number_cells(), 0.0)));
  }

  // Check for adjoint calculation
  if (d_input->check("adjoint"))
    d_adjoint = 0 != d_input->get<int>("adjoint");
}

//---------------------------------------------------------------------------//
void State::clear()
{
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      d_moments[g][i] = 0.0;
      if (d_store_angular_flux)
      {
        for (size_t a = 0; a < d_quadrature->number_angles(); ++a)
        {
          d_angular_flux[g][a][i] = 0.0;
        }
      }
    }
  }
}

//---------------------------------------------------------------------------//
void State::scale(const double f)
{
  for (size_t g = 0; g < d_number_groups; ++g)
  {
    for (size_t i = 0; i < d_mesh->number_cells(); ++i)
    {
      d_moments[g][i] *= f;
      if (d_store_angular_flux)
      {
        for (size_t a = 0; a < d_quadrature->number_angles(); ++a)
        {
          d_angular_flux[g][a][i] *= f;
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
  for (size_t g = 0; g < d_moments.size() + 1; g++)
    printf("--------------");
  printf("\n");
  printf("Scalar Flux Moments\n");
  printf("cell \\ g");
  for (size_t g = 0; g < d_number_groups; g++)
    printf(" %12i ", g);
  printf("\n");
  for (size_t g = 0; g < d_moments.size() + 1; g++)
    printf("--------------");
  printf("\n");
  for (size_t i = 0; i < d_mesh->number_cells(); i++)
  {
    printf("%10i", i);
    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf(" %12.5e ", d_moments[g][i]);
    }
    printf("\n");
  }

  if (d_store_angular_flux)
  {
    printf("\n");
    for (size_t a = 0; a < d_angular_flux[0].size() + 1; a++)
      printf("--------------");
    printf("\n");
    printf("Discrete Angular Flux\n");
    for (size_t a = 0; a < d_angular_flux[0].size() + 1; a++)
      printf("--------------");
    printf("\n");

    for (size_t g = 0; g < d_number_groups; g++)
    {
      printf("group %4i \n", g);
      printf("cell \\ a");
      for (size_t a = 0; a < d_angular_flux[g].size(); a++)
        printf(" %12i ", a);
      printf("\n");
      for (size_t a = 0; a < d_angular_flux[0].size() + 1; a++)
        printf("--------------");
      printf("\n");
      for (size_t i = 0; i < d_mesh->number_cells(); i++)
      {
        printf("%10i", i);
        for (size_t a = 0; a < d_angular_flux[g].size(); a++)
        {
          printf(" %12.5e ", d_angular_flux[g][a][i]);
        }
        printf("\n");
      }
      printf("\n");
    }
  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of State.cc
//---------------------------------------------------------------------------//
