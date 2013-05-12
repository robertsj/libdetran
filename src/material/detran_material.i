//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_material.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran materials.
 */
//---------------------------------------------------------------------------//

%module(package="detran") material
%{
#include <stddef.h>
#include "material/detran_material.hh"
%}

%feature("autodoc", "3");

// Hide templates from SWIG
%inline
{
#define MATERIAL_EXPORT
#define MATERIAL_TEMPLATE_EXPORT(...)
#define MATERIAL_INSTANTIATE_EXPORT(...)
}

%import "utilities/detran_utilities.i"

%include "Material.hh"

%template(MaterialSP) detran_utilities::SP<detran_material::Material>;

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//