//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_materials.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran kinetics.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"
%include "detran_material.i"

// Sources
%include "TimeDependentExternalSource.hh"
%include "LinearExternalSource.hh"
%include "PulsedExternalSource.hh"
%template(TDSourceSP) detran_utilities::SP<detran::TimeDependentExternalSource>;

%template(vec_source) std::vector<detran_utilities::SP<detran_external_source::ExternalSource> >;

// Materials
%include "KineticsMaterial.hh"
%include "TimeDependentMaterial.hh"
%include "LinearMaterial.hh"
%template(KineticsMaterialSP)   detran_utilities::SP<detran::KineticsMaterial>;
%template(TDMaterialSP)         detran_utilities::SP<detran::TimeDependentMaterial>;

%template(vec_material)         std::vector<detran_utilities::SP<detran::KineticsMaterial> >;

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//