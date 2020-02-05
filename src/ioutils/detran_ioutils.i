//----------------------------------*-C++-*-------------------------------0---//
/**
 *  @file   detran_ioutils.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran IO utilities.
 */
//----------------------------------------------------------------------------//

%module(directors="1", allprotected="1", package="detran") ioutils
%{
#include <stddef.h>
#include "ioutils/IO_HDF5.hh"
#include "ioutils/SiloOutput.hh"
#include "ioutils/PPMPlotter.hh"
#include "ioutils/PPMOutput.hh"
#include "ioutils/ColorMap.hh"
%}

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define IOUTILS_EXPORT
#define IOUTILS_TEMPLATE_EXPORT(...)
#define IOUTILS_INSTANTIATE_EXPORT(...)
}

%import "detran_transport.i"

%include "detran_config.hh"

#ifdef DETRAN_ENABLE_HDF5
%include "IO_HDF5.hh"
#endif

#ifdef DETRAN_ENABLE_SILO
%include "SiloOutput.hh"
#endif

%include "PPMPlotter.hh"
%include "PPMOutput.hh"
%include "ColorMap.hh"

//----------------------------------------------------------------------------//
//              end of detran_ioutils.i
//----------------------------------------------------------------------------//
