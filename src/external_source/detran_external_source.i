//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_external_source.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran external source.
 */
//---------------------------------------------------------------------------//

%include "ExternalSource.hh"
%include "ConstantSource.hh"
%include "DiscreteSource.hh"
%include "IsotropicSource.hh"

%template(ExternalSourceSP)  detran_utilities::SP<detran_external_source::ExternalSource>;
%template(ConstantSourceSP)  detran_utilities::SP<detran_external_source::ConstantSource>;
%template(DiscreteSourceSP)  detran_utilities::SP<detran_external_source::DiscreteSource>;
%template(IsotropicSourceSP) detran_utilities::SP<detran_external_source::IsotropicSource>;

//---------------------------------------------------------------------------//
//              end of detran_external_source.i
//---------------------------------------------------------------------------//





