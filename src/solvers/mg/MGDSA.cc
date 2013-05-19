//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MGDSA.cc
 *  @brief MGDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGDSA.hh"
#include "callow/solver/LinearSolverCreator.hh"
#include "ioutils/StdOutUtils.hh"

namespace detran
{

//----------------------------------------------------------------------------//
MGDSA::MGDSA(SP_input         input,
             SP_material      material,
             SP_mesh          mesh,
             SP_scattersource source,
             size_t           cutoff,
             bool             include_fission)
  : Base(input, material, mesh, cutoff, "MG-DSA")
  , d_scattersource(source)
{
  Require(d_scattersource);

  // Compute the diffusion coefficients.
  // \todo Need a flag that says whether they are build or not.
  d_material->compute_diff_coef();

  // Check for inner preconditioner db
  SP_input db;
  if (d_input->check("outer_pc_db"))
    db = d_input->get<SP_input>("outer_pc_db");

  // Create the loss operator for this group
  d_operator = new Operator_T(d_input,
                              d_material,
                              d_mesh,
                              include_fission,
                              d_group_cutoff,
                              false, // adjoint
                              1.0);  // keff


  // Create the linear solver for this group.
  d_solver = callow::LinearSolverCreator::Create(db);

  // Set the operators for this group.  The database is used
  // to set the preconditioner parameters for the diffusion solves.
  d_solver->set_operators(d_operator, db);

  // DEBUG -- inherit from MatrixShell
  // set_size(d_operator->number_columns());
  // d_operator->compute_explicit("mgdiff.out");
  // compute_explicit("mgdsa.out");
  // d_operator->print_matlab("mg.out");
}

//----------------------------------------------------------------------------//
void MGDSA::apply(Vector &V_in, Vector &V_out)
{

//  V_out.copy(V_in);
//  return;

  // Currently, DSA is only used on the flux moments,
  // not on the boundaries.
  size_t size_moments = d_mesh->number_cells();
  V_out.set(0.0);
  // Copy input vector to a multigroup flux; only the Krylov block is used.
  State::vec_moments_type
    phi(d_number_groups, State::moments_type(size_moments, 0.0));
  for (int g = d_group_cutoff; g < d_number_groups; ++g)
    for (int i = 0; i < size_moments; ++i)
      phi[g][i] = V_in[(g - d_group_cutoff) * size_moments + i];

  // Create the total group source, and copy into temporary Z
  Vector Z(size_moments * d_number_active_groups, 0.0);
  for (int g = d_group_cutoff; g < d_number_groups; g++)
  {
    State::moments_type source(size_moments, 0.0);
    d_scattersource->build_total_group_source(g, d_group_cutoff, phi, source);
    //detran_ioutils::print_vec(source);
    for (int i = 0; i < size_moments; i++)
      Z[(g - d_group_cutoff) * size_moments + i] = source[i];
  }
//  Z.print_matlab("Z.out");
  //-------------------------------------------------------------------------//
  // OPERATE: V2 <-- inv(C)*V1 = inv(C)*S*V0.
  //-------------------------------------------------------------------------//
  Vector Z_out(size_moments * d_number_active_groups, 0.0);
  d_solver->solve(Z, Z_out);
//  Z_out.print_matlab("Z2.out");
//  THROW("lala");
  //-------------------------------------------------------------------------//
  // OPERATE: V3 <-- V0 + V2 = V0 + inv(C)*S*V0 = (I + inv(C)*S)*V0
  //-------------------------------------------------------------------------//

  // \todo Implement an interface for swapping internal pointers a la PETSc
  V_out.copy(V_in);
  for (int i = 0; i < size_moments * d_number_active_groups; ++i)
    V_out[i] += Z_out[i];

}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGDSA.cc
//----------------------------------------------------------------------------//
