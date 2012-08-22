//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   PreconditionerMG.cc
 * \brief  PreconditionerMG 
 * \author Jeremy Roberts
 * \date   Aug 6, 2012
 */
//---------------------------------------------------------------------------//

#include "PreconditionerMG.hh"

namespace detran
{

PreconditionerMG::PreconditionerMG(SP_input input,
                                   SP_material material,
                                   SP_mesh mesh,
                                   SP_scattersource source,
                                   int moments_size_group,
                                   int boundary_size_group,
                                   int upscatter_cutoff,
                                   PC mg_pc)
  : d_input(input)
  , d_material(material)
  , d_mesh(mesh)
  , d_scattersource(source)
  , d_moments_size_group(moments_size_group)
  , d_boundary_size_group(boundary_size_group)
  , d_upscatter_cutoff(upscatter_cutoff)
  , d_tolerance(1.0e-5)
{
  // Preconditions
  Require(d_input);
  Require(d_material);
  Require(d_mesh);
  Require(d_scattersource);

  if (boundary_size_group)
    THROW("Reflective conditions for MG Preconditioner not yet supported");

  PetscErrorCode ierr;

  // Compute the diffusion coefficients. \todo There might be a better approach
  d_material->compute_diff_coef();

  // Set tolerance
  if (d_input->check("outer_pc_tolerance"))
    d_tolerance = input->get<double>("outer_pc_tolerance");

  // PC type
  std::string pc_type = "ilu";
  if (d_input->check("outer_pc_type"))
    pc_type = input->get<std::string>("outer_pc_type");

  // ILU levels
  int levels = 0;
  if (d_input->check("outer_pc_ilu_levels"))
    levels = input->get<int>("outer_pc_ilu_levels");

  // Create the diffusion operator
  d_lossoperator = new detran_diffusion::
    LossOperator(d_input, d_material, d_mesh);

  // Create the linear solver for this group.
  ierr = KSPCreate(PETSC_COMM_SELF, &d_solver);

  // Set the operators for this group.
  ierr = KSPSetOperators(d_solver,
                         d_lossoperator->get_operator(),
                         d_lossoperator->get_operator(),
                         SAME_PRECONDITIONER);

  // Set tolerances.
  ierr = KSPSetTolerances(d_solver,
                          d_tolerance,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT,
                          PETSC_DEFAULT);

  // Set the default preconditioner, ILU with a factor level of 2.
  // This has been in my experience a fair setting.
  PC diffusion_pc;
  KSPGetPC(d_solver, &diffusion_pc);

  PCSetType(diffusion_pc, pc_type.c_str());
  PCFactorSetLevels(diffusion_pc, levels);

  // Allow for command line flags.
  // \todo We need a way to separate the flags for these solvers and the
  // main transport solve.  For now, keep this as the default.
  ierr = KSPSetFromOptions(d_solver);

  // Set the PC context
  ierr = PCShellSetContext(mg_pc, this);
  Insist(!ierr, "Error setting within-group preconditioner context.");

  // Set the PC operator
  ierr = PCShellSetApply(mg_pc, apply_inv_P_MG);
  Insist(!ierr, "Error setting within-group preconditioner operator.");

  // Set the PC name for good measure
  ierr = PCShellSetName(mg_pc, "Multigroup Diffusion Preconditioner");
  Insist(!ierr, "Error within-group preconditioner name.");

  // Create the Vec's for use within the solver.  These are
  // sized to the solver, since the ones passed for applications
  // of the PC might be larger than the scalar flux vector due
  // to reflecting boundaries.  Here, only the fluxes are preconditioned.
  // Create the corresponding vectors.
  ierr = VecCreateSeq(PETSC_COMM_SELF,
                      d_material->number_groups() * d_moments_size_group,
                      &d_x);
  ierr = VecCreateSeq(PETSC_COMM_SELF,
                      d_material->number_groups() * d_moments_size_group,
                      &d_y);

}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PreconditionerMG.cc
//---------------------------------------------------------------------------//
