//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_materials.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran materials.
 */
//---------------------------------------------------------------------------//

%module detran_materials
%{
#include "SP.hh"
#include "Material.hh"
%}

// Load the standard library interfaces
%include std_vector.i

%include std_vec_typemap.i

%apply (std::vector<int> INPUTVECTOR)
       {(std::vector<int> value)}
%apply (std::vector<double> INPUTVECTOR)
       {(std::vector<double> value)}
       
%include "Definitions.hh"
%include "SP.hh"
%include "Material.hh"


%template(MaterialSP) detran_utils::SP<detran::Material>;
