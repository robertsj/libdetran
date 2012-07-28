//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_diffusion.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran diffusion.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"

#ifdef DETRAN_ENABLE_SLEPC
%include "DiffusionEigensolver.hh"
#endif

%include "LossOperator.hh"
%include "GainOperator.hh"

%template(LossOperatorSP) detran::SP<detran_diffusion::LossOperator>;
%template(GainOperatorSP) detran::SP<detran_diffusion::GainOperator>;

//---------------------------------------------------------------------------//
//              end of detran_diffusion.i
//---------------------------------------------------------------------------//
