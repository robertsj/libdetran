//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerWG.cc
 * \brief  PreconditionerWG 
 * \author Jeremy Roberts
 * \date   Aug 5, 2012
 */
//---------------------------------------------------------------------------//

#include "PreconditionerWG.hh"

#include "OneGroupLossOperator.hh"

namespace detran
{

PreconditionerWG::PreconditionerWG(SP_input input,
                                   SP_material material,
                                   SP_mesh mesh,
                                   SP_scattersource source,
                                   PC wg_pc)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_scattersource(source)
  , d_group(0)
  , d_solver(d_material->number_groups())
  , d_lossoperator(d_material->number_groups())
  , d_tolerance(1.0e-5)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_scattersource);

  PetscErrorCode ierr;

  // Compute the diffusion coefficients. \todo There might be a better approach
  d_material->compute_diff_coef();

  // Set tolerance
  if (d_input->check("inner_pc_tolerance"))
    d_tolerance = input->get<double>("inner_pc_tolerance");

  // PC type
  std::string pc_type = "ilu";
  if (d_input->check("inner_pc_type"))
    pc_type = input->get<std::string>("inner_pc_type");

  // ILU levels
  int levels = 0;
  if (d_input->check("inner_pc_ilu_levels"))
    levels = input->get<int>("inner_pc_ilu_levels");

  // Create the group-wise diffusion operators and the
  // associated linear systems for applying the inverse.
  for (int g = 0; g < d_material->number_groups(); g++)
  {

    // Create the loss operator for this group
    d_lossoperator[g] = new detran_diffusion::
      OneGroupLossOperator(d_input, d_material, d_mesh, g);

    // Create the linear solver for this group.
    ierr = KSPCreate(PETSC_COMM_SELF, &d_solver[g]);

    // Set the operators for this group.
    ierr = KSPSetOperators(d_solver[g],
                           d_lossoperator[g]->get_operator(),
                           d_lossoperator[g]->get_operator(),
                           SAME_PRECONDITIONER);

    // Set tolerances.
    ierr = KSPSetTolerances(d_solver[g],
                            d_tolerance,
                            PETSC_DEFAULT,
                            PETSC_DEFAULT,
                            PETSC_DEFAULT);

    // Set the default precondition, ILU with a factor level of 2.
    // This has been in my experience a fair setting.
    PC diffusion_pc;
    KSPGetPC(d_solver[g], &diffusion_pc);

    PCSetType(diffusion_pc, pc_type.c_str());
    PCFactorSetLevels(diffusion_pc, levels);

    // Allow for command line flags.
    // \todo We need a way to separate the flags for these solvers and the
    // main transport solve.  For now, keep this as the default.
    ierr = KSPSetFromOptions(d_solver[g]);
  }

  // Set the PC context
  ierr = PCShellSetContext(wg_pc, this);
  Insist(!ierr, "Error setting within-group preconditioner context.");

  // Set the PC operator
  ierr = PCShellSetApply(wg_pc, apply_inv_P);
  Insist(!ierr, "Error setting within-group preconditioner operator.");

  // Set the PC name for good measure
  ierr = PCShellSetName(wg_pc, "DSA Preconditioner");
  Insist(!ierr, "Error within-group preconditioner name.");

  // Create the Vec's for use within the solver.  These are
  // sized to the solver, since the ones passed for applications
  // of the PC might be larger than the scalar flux vector due
  // to reflecting boundaries.  Here, only the fluxes are preconditioned.
  // Create the corresponding vectors.
  ierr = VecCreateSeq(PETSC_COMM_SELF,
                      d_mesh->number_cells(),
                      &d_x);
  ierr = VecCreateSeq(PETSC_COMM_SELF,
                      d_mesh->number_cells(),
                      &d_y);

}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PreconditionerWG.cc
//---------------------------------------------------------------------------//
