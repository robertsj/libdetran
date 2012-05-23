//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   detran_utilities.i
 * \author Jeremy Roberts
 * \brief  Python interface for detran utilities.
 */
//---------------------------------------------------------------------------//

// STL
%include std_map.i
%include std_string.i
%include std_vector.i

// SP
%include "SP.hh"

// Vectors
%include "Definitions.hh"
namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

// Dimension
%include "Traits.hh"

// Input
%include "InputDB.hh"

// SP template of input database
%template(InputDBSP) detran::SP<detran::InputDB>;

// Template getters and setters
%template(get_int)      detran::InputDB::get<int>;
%template(get_dbl)      detran::InputDB::get<double>;
%template(get_vec_int)  detran::InputDB::get<detran::vec_int>;
%template(get_vec_dbl)  detran::InputDB::get<detran::vec_dbl>;
%template(get_str)      detran::InputDB::get<std::string>;
%template(put_int)      detran::InputDB::put<int>;
%template(put_dbl)      detran::InputDB::put<double>;
%template(put_vec_int)  detran::InputDB::put<detran::vec_int>;
%template(put_vec_dbl)  detran::InputDB::put<detran::vec_dbl>;
%template(put_str)      detran::InputDB::put<std::string>;

//---------------------------------------------------------------------------//
//              end of detran_utilities.i
//---------------------------------------------------------------------------//
