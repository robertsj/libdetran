//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGCMDSA.cc
 *  @brief MGCMDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "MGCMDSA.hh"

namespace detran
{

//----------------------------------------------------------------------------//
MGCMDSA::MGCMDSA(SP_input         input,
                 SP_material      material,
                 SP_mesh          mesh,
                 SP_scattersource source,
                 size_t           cutoff,
                 bool             include_fission)
  : Base(input, material, mesh, cutoff, "MGCMDSA")
{
  Require(d_scattersource);

  // Compute the diffusion coefficients.
  d_material->compute_diff_coef();

  // Check for inner preconditioner db
  SP_input db;
  if (d_input->check("outer_pc_db"))
    db = d_input->get<SP_input>("outer_pc_db");

  // Create the coarse mesh
  d_coarsemesher = new CoarseMesh(d_mesh, 1);

  // Homogenize the material
  Homogenize H(d_material);
  SP_state state(new State(d_input, d_mesh));
  SP_material coarse_material = H.homogenize(state, d_mesh, "COARSEMESH");

  // Create the loss operator for this group
  d_operator = new Operator_T(d_input,
                              coarse_material,
                              d_coarsemesher->get_coarse_mesh(),
                              include_fission,
                              d_group_cutoff,
                              false, // adjoint
                              1.0);  // keff

}

//----------------------------------------------------------------------------//
void MGCMDSA::apply(Vector &b, Vector &x)
{

}

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGCMDSA.cc
//----------------------------------------------------------------------------//
