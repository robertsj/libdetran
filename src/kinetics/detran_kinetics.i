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


// Downcasts and Upcasts for generic routines
%inline
{
  
  // Kinetics -> Base
  detran_utilities::SP<detran_material::Material> 
  downcast(detran_utilities::SP<detran::KineticsMaterial>* p)
  {
    return detran_utilities::SP<detran_material::Material>(*p);
  } 
  
  // TD -> Base
  detran_utilities::SP<detran_material::Material> 
  downcast(detran_utilities::SP<detran::TimeDependentMaterial>* p)
  {
    return detran_utilities::SP<detran_material::Material>(*p);
  } 
  
//  // Base -> Kinetics
//  detran_utilities::SP<detran::KineticsMaterial> 
//  downcast(detran_utilities::SP<detran_material::Material>* p)
//  {
//    return detran_utilities::SP<detran::KineticsMaterial>(*p);
//  } 
  
}

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//