//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGTransportSolver.cc
 *  @brief MGTransportSolver member definitions
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "detran_config.hh"
#include "MGTransportSolver.hh"
#include "WGSolverSI.hh"
#include "WGSolverGMRES.hh"

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
MGTransportSolver<D>::MGTransportSolver(SP_state                  state,
                                        SP_material               material,
                                        SP_boundary               boundary,
                                        const vec_externalsource &q_e,
                                        SP_fissionsource          q_f,
                                        bool                      multiply)
  : Base(state, material, boundary, q_e, q_f, multiply)
{

  // Get the quadrature from the state
  d_quadrature = d_state->get_quadrature();
  Ensure(d_quadrature);

  // Get the inner solver type and create.
  std::string wg_solver = "SI";
  if (d_input->check("inner_solver"))
  {
    wg_solver = d_input->template get<std::string>("inner_solver");
  }
  if (wg_solver == "SI")
  {
    d_wg_solver = new WGSolverSI<D>(d_state, d_material, d_quadrature,
                                    d_boundary, d_externalsources,
                                    d_fissionsource, d_multiply);
  }
  else if (wg_solver == "GMRES")
  {
    d_wg_solver = new WGSolverGMRES<D>(d_state, d_material, d_quadrature,
                                       d_boundary, d_externalsources,
                                       d_fissionsource, d_multiply);
  }
  else
  {
    THROW("Unsupported inner solver type selected: " + wg_solver);
  }

}

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//---------------------------------------------------------------------------//

template class MGTransportSolver<_1D>;
template class MGTransportSolver<_2D>;
template class MGTransportSolver<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of file MGTransportSolver.cc
//----------------------------------------------------------------------------//
