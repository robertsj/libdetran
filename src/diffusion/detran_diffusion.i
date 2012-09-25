//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_diffusion.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran diffusion.
 */
//---------------------------------------------------------------------------//

#ifdef DETRAN_ENABLE_PETSC
%include "BaseOperator.hh"
%include "LossOperator.hh"
%include "GainOperator.hh"
%template(LossOperatorSP) detran_utilities::SP<detran_diffusion::LossOperator>;
%template(GainOperatorSP) detran_utilities::SP<detran_diffusion::GainOperator>;
#endif

//---------------------------------------------------------------------------//
//              end of detran_diffusion.i
//---------------------------------------------------------------------------//
