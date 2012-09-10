//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_diffusion.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran diffusion.
 */
//---------------------------------------------------------------------------//

#ifdef DETRAN_ENABLE_SLEPC
%include "DiffusionEigensolver.hh"
#endif

%include "BaseOperator.hh"
%include "LossOperator.hh"
%include "GainOperator.hh"

%template(LossOperatorSP) detran_utilities::SP<detran_diffusion::LossOperator>;
%template(GainOperatorSP) detran_utilities::SP<detran_diffusion::GainOperator>;

//---------------------------------------------------------------------------//
//              end of detran_diffusion.i
//---------------------------------------------------------------------------//
