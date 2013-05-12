//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MGCMDSA.cc
 *  @brief MGCMDSA member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "MGCMDSA.hh"

namespace detran
{

//---------------------------------------------------------------------------//
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
  // \todo Need a flag that says whether they are build or not.
  d_material->compute_diff_coef();

  // Check for inner preconditioner db
  SP_input db;
  if (d_input->check("outer_pc_db"))
    db = d_input->get<SP_input>("outer_pc_db");



}

void MGCMDSA::apply(Vector &b, Vector &x)
{

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file MGCMDSA.cc
//---------------------------------------------------------------------------//
