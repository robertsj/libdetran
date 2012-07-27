//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_diffusion.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran diffusion.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"

%include "DiffusionEigensolver.hh"
%include "LossOperator.hh"
%include "GainOperator.hh"

//%include "detran_config.h"

//#ifdef DETRAN_ENABLE_PETSC

%template(LossOperatorSP) detran::SP<detran_diffusion::LossOperator>;
%template(GainOperatorSP) detran::SP<detran_diffusion::GainOperator>;

//#endif

//---------------------------------------------------------------------------//
//              end of detran_diffusion.i
//---------------------------------------------------------------------------//
