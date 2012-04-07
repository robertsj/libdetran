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
//%include "SourceIteration.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class SourceIteration
 * \brief 
 */
//===========================================================================//
template <class D>
class SourceIteration: public InnerIteration<D>
{
public:
  // SP Constructor
  static detran_utils::SP<SourceIteration<D> >
  Create(detran_utils::SP<detran_utils::InputDB>   input,
         detran_utils::SP<detran::State>           state,
         detran_utils::SP<detran::Mesh>            mesh,
         detran_utils::SP<detran::Material>        material,
         detran_utils::SP<detran::Quadrature>      quadrature,
         detran_utils::SP<detran::Boundary<D> >    boundary,
         detran_utils::SP<detran::ExternalSource>  q_e,
         detran_utils::SP<detran::FissionSource>   q_f);
         
  // Solve
  void solve(int g); 
};

} // end namespace detran


//%template(StateSP)  detran_utils::SP<detran::State>;

%template(InnerIteration2D)     detran::InnerIteration<detran::_2D>;
%template(InnerIteration2DSP)   detran_utils::SP<detran::InnerIteration<detran::_2D> >;

%template(SourceIteration2D)    detran::SourceIteration<detran::_2D>;
%template(SourceIteration2DSP)  detran_utils::SP<detran::SourceIteration<detran::_2D> >;


//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//





