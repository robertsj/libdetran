//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Preconditioner.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for callow preconditioners
 */
//---------------------------------------------------------------------------//

%include "detran_utilities.i"

%include "Preconditioner.hh"

%template(PreconditionerSP) detran_utilities::SP<callow::Preconditioner>;

