//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_transport.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran geometry.
 */
//---------------------------------------------------------------------------//

%module detran_transport
%{
#include "Definitions.hh"
#include "SP.hh"
#include "Quadrature.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "State.hh"
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
%include "State.hh"
%include "Quadrature.hh"
%include "Mesh.hh"

namespace detran
{



} // end namespace detran

%template(StateSP)  detran_utils::SP<detran::State>;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





