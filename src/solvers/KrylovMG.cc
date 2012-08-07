//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   KrylovMG.cc
 * \author robertsj
 * \date   Jun 19, 2012
 * \brief  GaussSeidel member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#include "detran_config.h"

#ifdef DETRAN_ENABLE_PETSC

// Detran
#include "KrylovMG.hh"

// System
#include <iostream>

namespace detran
{

// Constructor
template <class D>
KrylovMG<D>::KrylovMG(SP_input          input,
                      SP_state          state,
                      SP_mesh           mesh,
                      SP_material       material,
                      SP_quadrature     quadrature,
                      SP_boundary       boundary,
                      SP_externalsource q_e,
                      SP_fissionsource  q_f)
  : Base(input, state, mesh, material, quadrature, boundary, q_e, q_f)
  , d_moments_size(0)
  , d_moments_size_group(0)
  , d_boundary_size(0)
  , d_boundary_size_group(0)
{

  //-------------------------------------------------------------------------//
  // DETERMINE ENERGY GROUP BOUNDS
  //-------------------------------------------------------------------------//

  // Set the bounds for the downscatter GS portion and upscatter Krylov
  // portion.  The default is to use GS on the downscatter block and Krylov
  // on the upscatter block.
  d_upscatter_cutoff = material->upscatter_cutoff();
  if (input->check("outer_upscatter_cutoff"))
  {
    d_upscatter_cutoff = input->get<int>("outer_upscatter_cutoff");
    Insist((d_upscatter_cutoff >= 0) and
           (d_upscatter_cutoff <= material->upscatter_cutoff()),
           "Upscatter cutoff must be >= 0 and <= material upscatter cutoff");
  }
  d_upscatter_size = d_number_groups - d_upscatter_cutoff;

  //-------------------------------------------------------------------------//
  // SET UNKNOWN SIZES FOR KRYLOV SOLVE
  //-------------------------------------------------------------------------//

  d_moments_size_group = d_state->moments_size();
  d_moments_size = d_upscatter_size * d_moments_size_group;

  // Determine the sizes of any reflected boundary fluxes.
  for (int side = 0; side < 2*D::dimension; side++)
  {
    if (boundary->is_reflective(side))
    {
      // We only need to store only half of the unknowns.
      d_boundary_size_group += boundary->boundary_flux_size(side)/2;
    }
  }
  d_boundary_size = d_upscatter_size * d_boundary_size_group;

  //-------------------------------------------------------------------------//
  // SETUP SWEEPER FOR MULTIGROUP OPERATOR
  //-------------------------------------------------------------------------//

  d_sweeper = d_inner_solver->get_sweeper();
  d_sweepsource = d_inner_solver->get_sweepsource();

  //-------------------------------------------------------------------------//
  // SETUP PETSC SOLVER
  //-------------------------------------------------------------------------//

  // Total number of unknowns in block to be solved via Krylov
  int number_unknowns = d_moments_size + d_boundary_size;

  PetscErrorCode ierr;

  ierr = MatCreateShell(PETSC_COMM_WORLD,
                        number_unknowns,
                        number_unknowns,
                        number_unknowns,
                        number_unknowns,
                        this,
                        &d_operator);
  Insist(!ierr, "Error creating MR shell matrix.");

  // Create the corresponding vectors.
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unknowns,
                      &d_X);
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unknowns,
                      &d_B);

  VecSet(d_X, 0.0);
  VecSet(d_B, 0.0);

  // Set the operator.
  ierr = set_operation();
  Insist(!ierr, "Error setting matrix-vector operation.");

  // Create the KSP object.
  ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP object.");

  // Set the operator.
  KSPSetOperators(d_solver, d_operator, d_operator, SAME_NONZERO_PATTERN);

  // Preconditioner
  if (d_input->check("outer_use_pc"))
  {
    // \todo Why do I get a compile error if I use d_input instead?
    if (input->get<int>("outer_use_pc"))
    {
      std::cout << "Using preconditioned multigroup Krylov." << std::endl;
      d_use_pc = true;
    }
  }
  if (d_use_pc)
  {
    // Set up the shell preconditioner
    PC pc;
    KSPGetPC(d_solver, &pc);
    PCSetType(pc, PCSHELL);

    // Create the preconditioner.
    d_pc = new PreconditionerMG(d_input,
                                d_material,
                                d_mesh,
                                d_sweepsource->get_scatter_source(),
                                d_moments_size_group,
                                d_boundary_size_group,
                                d_upscatter_cutoff,
                                pc);

    ierr = KSPSetPCSide(d_solver, PC_LEFT);
  }

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver,
                          d_tolerance,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT);

  // Allow for command line flags.
  ierr = KSPSetFromOptions(d_solver);

}

template <class D>
PetscErrorCode KrylovMG<D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_MGTO_3D);
}
template <>
PetscErrorCode KrylovMG<_2D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_MGTO_2D);
}
template <>
PetscErrorCode KrylovMG<_1D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_MGTO_1D);
}


} // end namespace detran

#endif

//---------------------------------------------------------------------------//
//              end of KrylovMG.cc
//---------------------------------------------------------------------------//
