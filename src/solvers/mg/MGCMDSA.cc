//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CMMGDSA.cc
 *  @brief  CMMGDSA 
 *  @author Jeremy Roberts
 *  @date   Mar 26, 2013
 */
//---------------------------------------------------------------------------//

#include "CMMGDSA.hh"

namespace detran
{

//---------------------------------------------------------------------------//
CMMGDSA::CMMGDSA(SP_input         input,
                 SP_material      material,
                 SP_mesh          mesh,
                 SP_scattersource source,
                 size_t           cutoff,
                 bool             include_fission)
  : Base(input, material, mesh, cutoff, "CMMGDSA")
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

void CMMGDSA::apply(Vector &b, Vector &x)
{

}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of file CMMGDSA.cc
//---------------------------------------------------------------------------//
