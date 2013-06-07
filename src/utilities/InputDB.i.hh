//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   InputDB.i.hh
 *  @author Jeremy Roberts
 *  @date   Mar 25, 2012
 *  @brief  InputDB inline member definitions
 */
//---------------------------------------------------------------------------//


#ifndef detran_utilities_INPUTDB_I_HH_
#define detran_utilities_INPUTDB_I_HH_

#include <iostream>

namespace detran_utilities
{

//---------------------------------------------------------------------------//
inline bool InputDB::check(const std::string &key) const
{
  using std::cout;
  using std::endl;
  using std::string;

  // Check integers
  std::map<std::string, int>::const_iterator it1;
  if (d_data_int.size() > 0)
  {
    it1 = d_data_int.find(key);
    if (it1 != d_data_int.end())
      return true;
  }

  // Check doubles
  std::map<std::string, double>::const_iterator it2;
  if (d_data_dbl.size() > 0)
  {
    it2 = d_data_dbl.find(key);
    if (it2 != d_data_dbl.end())
      return true;
  }

  // Check integer vectors
  std::map<std::string, vec_int>::const_iterator it3;
  if (d_data_vec_int.size() > 0)
  {
    it3 = d_data_vec_int.find(key);
    if (it3 != d_data_vec_int.end())
      return true;
  }

  // Check double vectors
  std::map<std::string, vec_dbl>::const_iterator it4;
  if (d_data_vec_dbl.size() > 0)
  {
    it4 = d_data_vec_dbl.find(key);
    if (it4 != d_data_vec_dbl.end())
      return true;
  }

  // Check strings
  std::map<std::string, std::string>::const_iterator it5;
  if (d_data_str.size() > 0)
  {
    it5 = d_data_str.find(key);
    if (it5 != d_data_str.end())
      return true;
  }

  // Check db
  std::map<std::string, SP_input>::const_iterator it6;
  if (d_data_db.size() > 0)
  {
    it6 = d_data_db.find(key);
    if (it6 != d_data_db.end())
      return true;
  }

  return false;
}

//---------------------------------------------------------------------------//
// Get

template <>
inline int InputDB::get<int>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_int.find(key)->second;
}
template <>
inline double InputDB::get<double>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_dbl.find(key)->second;
}
template <>
inline vec_int InputDB::get<vec_int>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_vec_int.find(key)->second;
}
template <>
inline vec_dbl InputDB::get<vec_dbl>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_vec_dbl.find(key)->second;
}
template <>
inline std::string InputDB::get<std::string>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_str.find(key)->second;
}
template <>
inline InputDB::SP_input InputDB::get<InputDB::SP_input>(const std::string &key) const
{
  Requirev(check(key), key);
  return d_data_db.find(key)->second;
}

//---------------------------------------------------------------------------//
// Put

template <>
inline void InputDB::put(const std::string &key, const int value)
{
  d_data_int[key] = value;
  Ensurev(check(key), key);
}
template <>
inline void InputDB::put(const std::string &key, const double value)
{
  d_data_dbl[key] = value;
  Ensurev(check(key), key);
}
template <>
inline void InputDB::put(const std::string &key, const vec_int value)
{
  d_data_vec_int[key] = value;
  Ensurev(check(key), key);
}
template <>
inline void InputDB::put(const std::string &key, const vec_dbl value)
{
  d_data_vec_dbl[key] = value;
  Ensurev(check(key), key);
}
template <>
inline void InputDB::put(const std::string &key, const std::string value)
{
  d_data_str[key] = value;
  Ensurev(check(key), key);
}
template <>
inline void InputDB::put(const std::string &key, const SP_input value)
{
  d_data_db[key] = value;
  Ensurev(check(key), key);
}

//---------------------------------------------------------------------------//
// Get map

template <class T>
inline const std::map<std::string, T>&
InputDB::get_map()
{
  return d_data_int;
}
template <>
inline const std::map<std::string, double>&
InputDB::get_map<double>()
{
  return d_data_dbl;
}
template <>
inline const std::map<std::string, std::string>&
InputDB::get_map<std::string>()
{
  return d_data_str;
}
template <>
inline const std::map<std::string, vec_int>&
InputDB::get_map<vec_int>()
{
  return d_data_vec_int;
}
template <>
inline const std::map<std::string, vec_dbl>&
InputDB::get_map<vec_dbl>()
{
  return d_data_vec_dbl;
}
template <>
inline const std::map<std::string, InputDB::SP_input>&
InputDB::get_map<InputDB::SP_input>()
{
  return d_data_db;
}

} // end namespace detran_utilities

#endif /* detran_utilities_INPUTDB_I_HH_ */

//---------------------------------------------------------------------------//
//              end of InputDB.i.hh
//---------------------------------------------------------------------------//
