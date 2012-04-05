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
#include "InnerIteration.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "SourceIteration.hh"
#include "State.hh"
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
%include "InnerIteration.hh"
%include "SourceIteration.hh"

namespace detran
{


} // end namespace detran


//%template(StateSP)  detran_utils::SP<detran::State>;

%template(SourceIteration1D)  detran_utils::SP<detran::SourceIteration<detran::_1D> >;
%template(SourceIteration2D)  detran_utils::SP<detran::SourceIteration<detran::_2D> >;
%template(SourceIteration3D)  detran_utils::SP<detran::SourceIteration<detran::_3D> >;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





