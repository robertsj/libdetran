//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   WGSolverGMRES.cc
 *  @author robertsj
 *  @date   Apr 4, 2012
 *  @brief  WGSolverGMRES class definition.
 */
//---------------------------------------------------------------------------//

#include "WGSolverGMRES.hh"
#include "PC_DSA.hh"
#include "callow/solver/LinearSolverCreator.hh"

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
 , d_update_boundary_flux(false)
{
  Require(d_input);

  // Create operator
  d_operator = new WGTransportOperator<D>(d_state,
                                          d_boundary,
                                          d_sweeper,
                                          d_sweepsource);

  // Create temporary unknown and right hand size vectors
  d_x = new callow::Vector(d_operator->number_rows(), 0.0);
  d_b = new callow::Vector(d_operator->number_rows(), 0.0);

  // Get callow solver parameter database
  SP_input db;
  if (d_input->check("inner_solver_db"))
  {
    db = d_input->template get<SP_input>("inner_solver_db");
  }

  d_solver = callow::LinearSolverCreator::Create(db);
  Assert(d_solver);

  // Set the transport operator.  Note, no second db argument
  // is given, since that is for setting PC's.  We do that
  // explicitly below.
  d_solver->set_operators(d_operator);

  // Check for boundary update.  In GMRES, this is not the
  // default, since it implies extra sweeps that may not be needed.
  if (d_input->check("compute_boundary_flux"))
  {
    if (d_input->template get<int>("compute_boundary_flux"))
    {
      d_update_boundary_flux =
        d_input->template get<int>("compute_boundary_flux");
    }
  }
 // d_operator->compute_explicit("WGTO.out");
  //--------------------------------------------------------------------------//
  // PRECONDITIONER
  //--------------------------------------------------------------------------//

  /*
   *  Eventual Options:
   *    0. none
   *    1. fine mesh diffusion (= DSA)
   *    2. coarse mesh diffusion (= CMDSA)
   *    3. coarse mesh transport (= ~S2 on coarse grid)
   */

  std::string pc_type = "none";
  if (d_input->check("inner_pc_type"))
  {
    pc_type = d_input->template get<std::string>("inner_pc_type");
  }
  //std::cout << "Using WG-GMRES with PC-" << pc_type;

  size_t pc_side = callow::LinearSolver::LEFT;
  if (d_input->check("inner_pc_side"))
  {
    pc_side = d_input->template get<int>("inner_pc_side");
  }

  if (pc_type == "DSA")
  {
    d_pc = new PC_DSA(d_input,
                      d_material,
                      d_mesh,
                      d_sweepsource->get_scatter_source());
  }

  if (d_pc)
  {
    d_solver->set_preconditioner(d_pc, pc_side);
  }

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class WGSolverGMRES<_1D>;
template class WGSolverGMRES<_2D>;
template class WGSolverGMRES<_3D>;

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of WGSolverGMRES.cc
//---------------------------------------------------------------------------//
