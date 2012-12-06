//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   detran_material.i
 *  @author Jeremy Roberts
 *  @brief  Python interface for detran materials.
 */
//---------------------------------------------------------------------------//

%module(package="pydetran") material
%{
#include <stddef.h>
#include "material/detran_material.hh"
%}

%import "utilities/detran_utilities.i"

%include "Material.hh"

%template(MaterialSP) detran_utilities::SP<detran_material::Material>;

//---------------------------------------------------------------------------//
//              end of detran_material.i
//---------------------------------------------------------------------------//