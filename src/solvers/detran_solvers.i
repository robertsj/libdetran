//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_solvers.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%module detran_solvers
%{
// Detran
#include "Boundary.hh"
#include "ExternalSource.hh"
#include "FissionSource.hh"
#include "InnerIteration.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "SourceIteration.hh"
#include "State.hh"
// Utilities
#include "Definitions.hh"
#include "InputDB.hh"
#include "SP.hh"
%}

// Load the standard library interfaces.
%include std_vector.i
%include std_string.i


%include "SP.hh"
%include "InnerIteration.hh"
%include "SourceIteration.hh"

namespace detran
{


} // end namespace detran


//%template(StateSP)  detran_utils::SP<detran::State>;

//%template(InnerIteration2D)     detran::InnerIteration<detran::_2D>;
//%template(InnerIteration2DSP)   detran_utils::SP<detran::InnerIteration<detran::_2D> >;
%template(SourceIteration1DSP)  detran_utils::SP<detran::SourceIteration<detran::_1D> >;
%template(SourceIteration2D)    detran::SourceIteration<detran::_2D>;
%template(SourceIteration2DSP)  detran_utils::SP<detran::SourceIteration<detran::_2D> >;
%template(SourceIteration3DSP)  detran_utils::SP<detran::SourceIteration<detran::_3D> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





