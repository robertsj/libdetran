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
  %template(Vec_int)  vector<int>;
  %template(Vec2_int) vector<vector<int> >;
  %template(Vec3_int) vector<vector<vector<int> > >;
  %template(Vec_dbl)  vector<double>;
  %template(Vec2_dbl) vector<vector<double> >;
  %template(Vec3_dbl) vector<vector<vector<double> > >;
}

// Dimension
%include "Traits.hh"

// Input
%include "InputDB.hh"

// Math utilities
%include "MathUtilities.hh"

// 2-D point
%include "Point.hh"

// SP template of input database
%template(InputDBSP) detran_utilities::SP<detran_utilities::InputDB>;

// Template getters and setters
%template(get_int)      detran_utilities::InputDB::get<int>;
%template(get_dbl)      detran_utilities::InputDB::get<double>;
%template(get_vec_int)  detran_utilities::InputDB::get<detran_utilities::vec_int>;
%template(get_vec_dbl)  detran_utilities::InputDB::get<detran_utilities::vec_dbl>;
%template(get_str)      detran_utilities::InputDB::get<std::string>;
%template(put_int)      detran_utilities::InputDB::put<int>;
%template(put_dbl)      detran_utilities::InputDB::put<double>;
%template(put_vec_int)  detran_utilities::InputDB::put<detran_utilities::vec_int>;
%template(put_vec_dbl)  detran_utilities::InputDB::put<detran_utilities::vec_dbl>;
%template(put_str)      detran_utilities::InputDB::put<std::string>;

//---------------------------------------------------------------------------//
//              end of detran_utilities.i
//---------------------------------------------------------------------------//
