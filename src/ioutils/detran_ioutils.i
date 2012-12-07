//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_ioutils.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran IO utilities.
 */
//---------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") ioutils
%{
#include <stddef.h>
#include "ioutils/IO_HDF5.hh"
#include "ioutils/SiloOutput.hh"
%}

%import "detran_transport.i"

%include "detran_config.hh"

#ifdef DETRAN_ENABLE_HDF5
%include "IO_HDF5.hh"
#endif

#ifdef DETRAN_ENABLE_SILO
%include "SiloOutput.hh"
#endif

//---------------------------------------------------------------------------//
//              end of detran_ioutils.i
//---------------------------------------------------------------------------//