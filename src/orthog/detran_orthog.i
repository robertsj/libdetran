//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_orthog.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran orthogonal basis sets.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") orthog
%{
#include <stddef.h>
#include "orthog/OrthogonalBasisParameters.hh"
#include "orthog/OrthogonalBasis.hh"
%}

// Hide templates from SWIG
%inline
{
#define ORTHOG_EXPORT
#define ORTHOG_TEMPLATE_EXPORT(...)
#define ORTHOG_INSTANTIATE_EXPORT(...)
}

%import "utilities/detran_utilities.i"
%import "callow/detran_callow.i"

// Base angle classes and utilities
%include "OrthogonalBasisParameters.hh"
%template(OrthogonalBasisParamsSP) detran_utilities::SP<detran_orthog::OrthogonalBasisParameters>;

%include "OrthogonalBasis.hh"
%template(OrthogonalBasisSP) detran_utilities::SP<detran_orthog::OrthogonalBasis>;

//---------------------------------------------------------------------------//
//              end of detran_orthog.i
//---------------------------------------------------------------------------//
