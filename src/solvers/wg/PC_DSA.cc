//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  PC_DSA.cc
 *  @brief PC_DSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "PC_DSA.hh"
#include "callow/solver/LinearSolverCreator.hh"

namespace detran
{

//---------------------------------------------------------------------------//
PC_DSA::PC_DSA(SP_input input,
               SP_material material,
               SP_mesh mesh,
               SP_scattersource source)
  : Base(input, material, mesh, "WG-DSA")
  , d_scattersource(source)
{
  // Preconditions
  Require(d_scattersource);

  // Compute the diffusion coefficients.
  // \todo Need a flag that says whether they are build or not.
  d_material->compute_diff_coef();

  // Check for inner preconditioner db
  SP_input db;
  if (d_input->check("inner_pc_db"))
    db = d_input->get<SP_input>("inner_pc_db");

  // Create the group-wise diffusion operators and the
  // associated linear systems for applying the inverse.
  for (int g = 0; g < d_number_groups; g++)
  {

    // Create the loss operator for this group
    d_operator[g] = new Operator_T(d_input, d_material, d_mesh, g);

    // Create the linear solver for this group.
    d_solver[g] = callow::LinearSolverCreator::Create(db);

    // Set the operators for this group.  The database is used
    // to set the preconditioner parameters for the diffusion solves.
    d_solver[g]->set_operators(d_operator[g], db);
  }

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file PC_DSA.cc
//---------------------------------------------------------------------------//
