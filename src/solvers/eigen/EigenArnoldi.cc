//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EigenArnoldi.cc
 *  @brief EigenArnoldi class definition.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "EigenArnoldi.hh"
#include "callow/solver/EigenSolverCreator.hh"

#include <iostream>

namespace detran
{

//----------------------------------------------------------------------------//
template <class D>
EigenArnoldi<D>::EigenArnoldi(SP_mg_solver mg_solver)
  : Base(mg_solver)
{

  // Create operator
  d_operator = new Operator_T(mg_solver);

  // Create vectors
  d_x  = new callow::Vector(d_operator->number_rows(), 0.0);
  d_x0 = new callow::Vector(d_operator->number_rows(), 0.0);

  // Get callow solver parameter database
  SP_input db;
  if (d_input->check("eigen_solver_db"))
  {
    db = d_input->template get<SP_input>("eigen_solver_db");
  }
  d_eigensolver = callow::EigenSolverCreator::Create(db);
  Assert(d_eigensolver);

  // Set the transport operator.
  d_eigensolver->set_operators(d_operator);
}

//----------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS
//----------------------------------------------------------------------------//

template class EigenArnoldi<_1D>;
template class EigenArnoldi<_2D>;
template class EigenArnoldi<_3D>;

} // end namespace detran

//----------------------------------------------------------------------------//
//              end of EigenArnoldi.cc
//----------------------------------------------------------------------------//
