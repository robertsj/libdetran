//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Solver.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow solvers
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"
%include "preconditioner/Preconditioner.i"

%include "LinearSolver.hh"
%include "EigenSolver.hh"
%include "LinearSolverCreator.hh"
%include "EigenSolverCreator.hh"

%template(LinearSolverSP) detran_utilities::SP<callow::LinearSolver>;
%template(EigenSolverSP)  detran_utilities::SP<callow::EigenSolver>;
