//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   WGSolverGMRES.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  WGSolverGMRES class definition.
 */
//---------------------------------------------------------------------------//

#include "detran_config.hh"

#ifdef DETRAN_ENABLE_PETSC

#include "WGSolverGMRES.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
WGSolverGMRES<D>::WGSolverGMRES(SP_state                  state,
                                SP_material               material,
                                SP_quadrature             quadrature,
                                SP_boundary               boundary,
                                const vec_externalsource &q_e,
                                SP_fissionsource          q_f,
                                bool                      multiply)
 : Base(state, material, quadrature, boundary, q_e, q_f, multiply)
 , d_moments_size(0)
 , d_boundary_size(0)
 , d_use_pc(false)
{

  // Error from PETSc calls.
  PetscErrorCode ierr;

  // Determine the sizes of the moments.
  d_moments_size = d_mesh->number_cells();

  // Determine the sizes of any reflected boundary fluxes.
  for (int side = 0; side < 2*D::dimension; side++)
  {
    if (boundary->is_reflective(side))
    {
      // We only need to store only half of the unknowns.
      d_boundary_size += boundary->boundary_flux_size(side)/2;
    }
  }

  // Total number of unknowns in the within-group solve.
  int number_unknowns = d_moments_size + d_boundary_size;

  ierr = MatCreateShell(PETSC_COMM_WORLD,
                        number_unknowns,
                        number_unknowns,
                        number_unknowns,
                        number_unknowns,
                        this,
                        &d_operator);
  Insist(!ierr, "Error creating transport operator shell matrix.");

  // Create the corresponding vectors.
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unknowns,
                      &d_X);
  ierr = VecCreateSeq(PETSC_COMM_WORLD,
                      number_unknowns,
                      &d_B);

  VecSet(d_X, 0.0);
  VecSet(d_B, 0.0);

  // Set the matrix-vector operator.
  ierr = set_operation();
  Insist(!ierr, "Error setting matrix-vector operation.");

  // Create the KSP object.
  ierr = KSPCreate(PETSC_COMM_WORLD, &d_solver);
  Insist(!ierr, "Error creating KSP object.");

  // Set the operator.
  KSPSetOperators(d_solver, d_operator, d_operator, SAME_NONZERO_PATTERN);

  //--------------------------------------------------------------------------//
  // PRECONDITIONER
  //--------------------------------------------------------------------------//

  /*
   *  Eventual Options:
   *    1. fine mesh diffusion (= DSA)
   *    2. coarse mesh diffusion (= CMDSA)
   *    3. coarse mesh transport (= ~S2 on coarse grid)
   */
  if (d_input->check("inner_use_pc"))
  {
    if (d_input->template get<int>("inner_use_pc"))
    {
      std::cout << "Using preconditioned GMRES." << std::endl;
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
    d_pc = new PreconditionerWG(d_input,
        d_material, d_mesh, d_sweepsource->get_scatter_source(), pc);

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
PetscErrorCode WGSolverGMRES<D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_3D);
}
template <>
PetscErrorCode WGSolverGMRES<_2D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_2D);
}
template <>
PetscErrorCode WGSolverGMRES<_1D>::set_operation()
{
  return MatShellSetOperation(d_operator,
                              MATOP_MULT,
                              (void(*)(void)) apply_WGTO_1D);
}

template class WGSolverGMRES<_1D>;
template class WGSolverGMRES<_2D>;
template class WGSolverGMRES<_3D>;

} // end namespace detran

#endif

//---------------------------------------------------------------------------//
//              end of WGSolverGMRES.cc
//---------------------------------------------------------------------------//
