//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_transport.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran transport.
 */
//---------------------------------------------------------------------------//

%module detran_transport
%{
// Detran
#include "Boundary.hh"
#include "FissionSource.hh"
#include "ExternalSource.hh"
#include "ConstantSource.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"
#include "Traits.hh"
// Utilities
#include "Definitions.hh"
#include "SP.hh"
%}

// Load the standard library interfaces.
%include std_vector.i
%include std_string.i

// Load the vector maps.  Note, I am a bit unhappy with
// how it's all used.  They work *if* I declare the class
// interface below.  Otherwise, just including e.g.
// Mesh2D doesn't allow the maps, since I'm using 
// typedefs on the input arguments.  There should be an
// easy way around this, but I'm not a SWIG pro.
%include std_vec_typemap.i

%include "SP.hh"

%include "Boundary.hh"
%include "FissionSource.hh"

// External/general sources
%include "ExternalSource.hh"
%include "ConstantSource.hh"
%include "State.hh"
%include "Traits.hh"
//%include "Mesh.hh"
namespace detran
{


} // end namespace detran

// Traits
//%template(_1D)

// Templates 

%template(StateSP)  detran_utils::SP<detran::State>;

%template(FissionSourceSP)  detran_utils::SP<detran::FissionSource>;

%template(ExternalSourceSP) detran_utils::SP<detran::ExternalSource>;
%template(ConstantSourceSP) detran_utils::SP<detran::ConstantSource>;

%template(Boundary1D)    detran::Boundary<detran::_1D>;
%template(Boundary1DSP)  detran_utils::SP<detran::Boundary<detran::_1D> >;
%template(Boundary2D)    detran::Boundary<detran::_2D>;
%template(Boundary2DSP)  detran_utils::SP<detran::Boundary<detran::_2D> >;
%template(Boundary3D)    detran::Boundary<detran::_3D>;
%template(Boundary3DSP)  detran_utils::SP<detran::Boundary<detran::_3D> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





