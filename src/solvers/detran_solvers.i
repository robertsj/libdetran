//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_solvers.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

// Note, only the basic, non-PETSc classes are exposed directly.
// The others can be accessed from within Execute.

%include "InnerIteration.hh"
%template(InnerIteration1D)     detran::InnerIteration<detran::_1D>;
%template(InnerIteration1DSP)   detran_utilities::SP<detran::InnerIteration<detran::_1D> >;
%template(InnerIteration2D)     detran::InnerIteration<detran::_2D>;
%template(InnerIteration2DSP)   detran_utilities::SP<detran::InnerIteration<detran::_2D> >;
%template(InnerIteration3D)     detran::InnerIteration<detran::_3D>;
%template(InnerIteration3DSP)   detran_utilities::SP<detran::InnerIteration<detran::_3D> >;

namespace detran
{

// This is the one class I couldn't easily wrap.  I narrow it down to it being
// derived and my subsequent use of the "using Base::something" syntax.  SWIG,
// for whatever reason, tries to make getters and setters for those base 
// variables, and the generated code fails due to typedef issues. 
template <class D>
class SourceIteration: public InnerIteration<D>
{
public:
  // Constructor
  SourceIteration(
      detran_utilities::SP<detran_utilities::InputDB>             input,
      detran_utilities::SP<detran::State>                           state,
      detran_utilities::SP<detran_geometry::Mesh>                   mesh,
      detran_utilities::SP<detran_material::Material>               material,
      detran_utilities::SP<detran_angle::Quadrature>                quadrature,
      detran_utilities::SP<detran::BoundaryBase<D> >                boundary,
      detran_utilities::SP<detran_external_source::ExternalSource>  q_e,
      detran_utilities::SP<detran::FissionSource>                   q_f);
  // SP Constructor
  static detran_utilities::SP<SourceIteration<D> >
  Create(
      detran_utilities::SP<detran_utilities::InputDB>             input,
      detran_utilities::SP<detran::State>                           state,
      detran_utilities::SP<detran_geometry::Mesh>                   mesh,
      detran_utilities::SP<detran_material::Material>               material,
      detran_utilities::SP<detran_angle::Quadrature>                quadrature,
      detran_utilities::SP<detran::BoundaryBase<D> >                boundary,
      detran_utilities::SP<detran_external_source::ExternalSource>  q_e,
      detran_utilities::SP<detran::FissionSource>                   q_f);
  // Solve the within group equation.
  void solve(const detran_utilities::size_t g);
};

} // end namespace detran

%template(SourceIteration1D)    detran::SourceIteration<detran::_1D>;
%template(SourceIteration1DSP)  detran_utilities::SP<detran::SourceIteration<detran::_1D> >;
%template(SourceIteration2D)    detran::SourceIteration<detran::_2D>;
%template(SourceIteration2DSP)  detran_utilities::SP<detran::SourceIteration<detran::_2D> >;


%include "SolverMG.hh"
%template(SolverMG1D)     detran::SolverMG<detran::_1D>;
%template(SolverMG2D)     detran::SolverMG<detran::_2D>;
%template(SolverMG3D)     detran::SolverMG<detran::_3D>;


%include "GaussSeidelMG.hh"
%template(GaussSeidelMG1D)        detran::GaussSeidelMG<detran::_1D>;
%template(GaussSeidelMG1DSP)      detran_utilities::SP<detran::GaussSeidelMG<detran::_1D> >;
%template(GaussSeidelMG2D)        detran::GaussSeidelMG<detran::_2D>;
%template(GaussSeidelMG2DSP)      detran_utilities::SP<detran::GaussSeidelMG<detran::_2D> >;


%include "Eigensolver.hh"
%template(Eigensolver1D)     detran::Eigensolver<detran::_1D>;
%template(Eigensolver2D)     detran::Eigensolver<detran::_2D>;
%template(Eigensolver3D)     detran::Eigensolver<detran::_3D>;

%include "PowerIteration.hh"
%template(PowerIteration1D)     detran::PowerIteration<detran::_1D>;
%template(PowerIteration1DSP)   detran_utilities::SP<detran::PowerIteration<detran::_1D> >;
%template(PowerIteration2D)     detran::PowerIteration<detran::_2D>;
%template(PowerIteration2DSP)   detran_utilities::SP<detran::PowerIteration<detran::_2D> >;
%template(PowerIteration3D)     detran::PowerIteration<detran::_3D>;
%template(PowerIteration3DSP)   detran_utilities::SP<detran::PowerIteration<detran::_3D> >;
//
//%include "DiffusionGainOperator.hh"
//%include "DiffusionLossOperator.hh"
//%include "DiffusionFixedSourceSolver.hh"
//%template(DiffusionFixed1D)     detran::DiffusionFixedSourceSolver<detran::_1D>;
//%template(DiffusionFixed1DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_1D> >;
//%template(DiffusionFixed2D)     detran::DiffusionFixedSourceSolver<detran::_2D>;
//%template(DiffusionFixed2DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_2D> >;
//%template(DiffusionFixed3D)     detran::DiffusionFixedSourceSolver<detran::_3D>;
//%template(DiffusionFixed3DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_3D> >;
//

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//
