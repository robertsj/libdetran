//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_materials.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran kinetics.
 */
//---------------------------------------------------------------------------//
   
%include "detran_utilities.i"
%include "detran_material.i"

// Include callback macros
%include callback.i
setCallbackMethod(1, PyTimeDependentMaterial, set_update_impl, (void* data), (), (), 1)


// Sources
%include "TimeDependentExternalSource.hh"
%include "LinearExternalSource.hh"
%include "PulsedExternalSource.hh"
%template(TDSourceSP) detran_utilities::SP<detran::TimeDependentExternalSource>;
%template(vec_source) std::vector<detran_utilities::SP<detran_external_source::ExternalSource> >;

// Materials
%include "KineticsMaterial.hh"
%include "TimeDependentMaterial.hh"
%include "PyTimeDependentMaterial.hh"
%include "LinearMaterial.hh"
%template(KineticsMaterialSP)   detran_utilities::SP<detran::KineticsMaterial>;
%template(TDMaterialSP)         detran_utilities::SP<detran::TimeDependentMaterial>;
%template(PYTDMaterialSP)       detran_utilities::SP<detran::PyTimeDependentMaterial>;
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

// Following construction of a Python material, set the 
// callback function
//%feature("pythonappend") detran::PyTimeDependentMaterial(const detran_utilities::size_t number_materials,
//                                                         const detran_utilities::size_t number_energy_groups,
//                                                         const detran_utilities::size_t number_precursor_groups,
//                                                         std::string  name) 
//%{
//   # the client *must* define update_impl
//   print "PYTHON setting implementation"
//   self.set_update_impl(self.update_impl)
//%}

// Set Python material update
//setCallbackMethod(1, PyTimeDependentMaterial, set_update_impl, (void* data), (), (), 1)

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//
