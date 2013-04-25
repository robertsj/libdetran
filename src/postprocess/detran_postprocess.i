//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_postprocess.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran postprocessing.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") postprocess
%{
#include <stddef.h>
#include "postprocess/ReactionRates.hh"
%}

// Hide templates from SWIG
%inline
{
#define POSTPROCESS_EXPORT
#define POSTPROCESS_TEMPLATE_EXPORT(...)
#define POSTPROCESS_INSTANTIATE_EXPORT(...)
}

%import "detran_transport.i"

%include "postprocess/ReactionRates.hh"

//---------------------------------------------------------------------------//
//              end of detran_postprocess.i
//---------------------------------------------------------------------------//