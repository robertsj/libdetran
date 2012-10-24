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
%include "Jacobi.hh"
%include "GaussSeidel.hh"
%include "GMRES.hh"
%include "Richardson.hh"

%include "EigenSolver.hh"
%include "PowerIteration.hh"

%template(LinearSolverSP) detran_utilities::SP<callow::LinearSolver>;
%template(EigenSolverSP)  detran_utilities::SP<callow::EigenSolver>;
