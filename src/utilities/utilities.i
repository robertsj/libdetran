//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   utilities.hh
 * \author Jeremy Roberts
 * \brief  Python interface for detran utilities.
 */
//---------------------------------------------------------------------------//

%module detran_utilities
%{
#include "Definitions.hh"
#include "DBC.hh"
#include "InputDB.hh"
%}

// Load the standard library interfaces
%include std_map.i
%include std_string.i
%include std_vector.i

%include std_vec_typemap.i
//%apply (int DIM1,  int* IN_ARRAY1)
//      {(int len1,  int* x1),
//       (int len2,  int* x2)}

%apply (std::vector<int> INPUTVECTOR)
       {(std::vector<int> value)}
%apply (std::vector<double> INPUTVECTOR)
       {(std::vector<double> value)}


%include "Definitions.hh"
%include "InputDB.hh"

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

namespace detran_utils
{

%template(get_int)      InputDB::get<int>;
%template(get_dbl)      InputDB::get<double>;
%template(get_vec_int)  InputDB::get<vec_int>;
%template(get_vec_dbl)  InputDB::get<vec_dbl>;
%template(get_str)      InputDB::get<std::string>;
%template(put_int)      InputDB::put<int>;
%template(put_dbl)      InputDB::put<double>;
%template(put_vec_int)  InputDB::put<vec_int>;
%template(put_vec_dbl)  InputDB::put<vec_dbl>;
%template(put_str)      InputDB::put<std::string>;

} // end namespace detran_utils






