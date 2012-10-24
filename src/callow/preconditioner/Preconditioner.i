//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Solver.i
 * \author Jeremy Roberts
 * \brief  Python interface for callow solvers
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Preconditioner.hh"
%include "PCILU0.hh"
%include "PCJacobi.hh"
%include "PCIdentity.hh"

%template(PreconditionerSP) detran_utilities::SP<callow::Preconditioner>;

