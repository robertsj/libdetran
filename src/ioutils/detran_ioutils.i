//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_ioutils.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran IO utilities.
 */
//---------------------------------------------------------------------------//


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