//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_solvers.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran solvers.
 */
//---------------------------------------------------------------------------//

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
