//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BaseOperator.cc
 * \brief  BaseOperator 
 * \author Jeremy Roberts
 * \date   Jul 25, 2012
 */
//---------------------------------------------------------------------------//

// Configuration
#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "BaseOperator.hh"

// Utilities
#include "DBC.hh"

#include <iostream>

namespace detran_diffusion
{

BaseOperator::BaseOperator(SP_input input,
                           SP_material material,
                           SP_mesh     mesh)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
{
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  d_dimension = d_mesh->dimension();
}

BaseOperator::~BaseOperator()
{
  // Clean up
  MatDestroy(&d_operator);
}

void BaseOperator::multiply(Vec v_in, Vec v_out)
{
  // Start timer
  d_timer.function_tic();

  // Preconditions
  //Require(v_in);
  //Require(v_out);

  // Error flag
  PetscErrorCode ierr;

  // Multiply
  ierr = MatMult(d_operator, v_in, v_out);
  //Assert(ierr);

  // End timer
  d_timer.function_toc("multiply");
}

void BaseOperator::display()
{
  // Error flag
  PetscErrorCode ierr;

  // Print to standard out.
  ierr = MatView(d_operator,  PETSC_VIEWER_STDOUT_SELF);
}

#endif

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file BaseOperator.cc
//---------------------------------------------------------------------------//
