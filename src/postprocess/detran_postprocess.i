//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_postprocess.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran postprocessing.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="pydetran") postprocess
%{
#include "postprocess/ReactionRates.hh"
%}

//---------------------------------------------------------------------------//
//              end of detran_postprocess.i
//---------------------------------------------------------------------------//