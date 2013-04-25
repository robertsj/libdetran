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
#include "orthog/OrthogonalBasis.hh"
#include "orthog/ContinuousOrthogonalBasis.hh"
#include "orthog/CLP.hh"
#include "orthog/DLP.hh"
#include "orthog/DCP.hh"
#include "orthog/DCT.hh"
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
%include "OrthogonalBasis.hh"
%include "ContinuousOrthogonalBasis.hh"
%include "CLP.hh"
%include "DLP.hh"
%include "DCP.hh"
%include "DCT.hh"

%template(OrthogonalBasisSP) detran_utilities::SP<detran_orthog::OrthogonalBasis>;

//---------------------------------------------------------------------------//
//              end of detran_orthog.i
//---------------------------------------------------------------------------//
