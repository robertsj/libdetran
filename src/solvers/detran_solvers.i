//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_solvers.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

%include "DiffusionGainOperator.hh"
%template(DiffusionGainOperatorSP)   detran_utilities::SP<detran::DiffusionGainOperator>;

%include "DiffusionLossOperator.hh"
%template(DiffusionLossOperatorSP)   detran_utilities::SP<detran::DiffusionLossOperator>;

%include "DiffusionFixedSourceSolver.hh"
%template(DiffusionFixed1D)     detran::DiffusionFixedSourceSolver<detran::_1D>;
%template(DiffusionFixed1DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_1D> >;
%template(DiffusionFixed2D)     detran::DiffusionFixedSourceSolver<detran::_2D>;
%template(DiffusionFixed2DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_2D> >;
%template(DiffusionFixed3D)     detran::DiffusionFixedSourceSolver<detran::_3D>;
%template(DiffusionFixed3DSP)   detran_utilities::SP<detran::DiffusionFixedSourceSolver<detran::_3D> >;

%include "DiffusionEigensolver.hh"
%template(DiffusionEigen1D)     detran::DiffusionEigensolver<detran::_1D>;
%template(DiffusionEigen1DSP)   detran_utilities::SP<detran::DiffusionEigensolver<detran::_1D> >;
%template(DiffusionEigen2D)     detran::DiffusionEigensolver<detran::_2D>;
%template(DiffusionEigen2DSP)   detran_utilities::SP<detran::DiffusionEigensolver<detran::_2D> >;
%template(DiffusionEigen3D)     detran::DiffusionEigensolver<detran::_3D>;
%template(DiffusionEigen3DSP)   detran_utilities::SP<detran::DiffusionEigensolver<detran::_3D> >;

%include "FixedSourceManager.hh"
%template(Fixed1D) detran::FixedSourceManager<detran::_1D>;
%template(Fixed2D) detran::FixedSourceManager<detran::_2D>;
%template(Fixed3D) detran::FixedSourceManager<detran::_3D>;

%include "EigenvalueManager.hh"
%template(Eigen1D) detran::EigenvalueManager<detran::_1D>;
%template(Eigen2D) detran::EigenvalueManager<detran::_2D>;
%template(Eigen3D) detran::EigenvalueManager<detran::_3D>;

//---------------------------------------------------------------------------//
//              end of detran_transport.i
//---------------------------------------------------------------------------//
