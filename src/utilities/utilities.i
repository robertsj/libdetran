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
#include "GenDB.hh"
%}

// Load the standard library interfaces
%include std_map.i
%include std_string.i
%include std_vector.i
// Load typemaps for vectors
%include std_vec_typemap.i
%apply (std::vector<int> INPUTVECTOR)
       {(std::vector<int> value)}
%apply (std::vector<double> INPUTVECTOR)
       {(std::vector<double> value)}

%include "Definitions.hh"

namespace std
{
  %template(vec_int) vector<int>;
  %template(vec_dbl) vector<double>;
}

namespace detran_utils
{

class GenDB 
{

public:

  /*
   *  \brief Constructor.
   */
  GenDB();

  /*
   *  \brief Return value of key.
   *  \param    key     Name of the parameter.
   *  \return           Check whether key is found.
   */
  template <class T>
  T get(const std::string &key) const;

  /*
   *  \brief Put a key and value in the database.
   *  \param    key         Name of the parameter.
   *  \param    value       Reference to which parameter value is assigned.
   *  \param    replace     Can we replace a current value?
   *  \return               Check whether key is found.
   */
  template <class T>
  void put(const std::string &key, const T value);

};


// Instantiate the important ones for Python use.
%template(get_int)      GenDB::get<int>;
%template(get_dbl)      GenDB::get<double>;
%template(get_vec_int)  GenDB::get<vec_int>;
%template(get_vec_dbl)  GenDB::get<vec_dbl>;
%template(get_str)      GenDB::get<std::string>;
%template(put_int)      GenDB::put<int>;
%template(put_dbl)      GenDB::put<double>;
%template(put_vec_int)  GenDB::put<vec_int>;
%template(put_vec_dbl)  GenDB::put<vec_dbl>;
%template(put_str)      GenDB::put<std::string>;

} // end namespace detran_utils






