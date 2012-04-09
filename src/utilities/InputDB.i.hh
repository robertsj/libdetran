//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputDB.i.hh
 * \author Jeremy Roberts
 * \date   Mar 25, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//


#ifndef INPUTDB_I_HH_
#define INPUTDB_I_HH_

#include <iostream>

namespace detran
{


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
  return false;
}

// Get

template <>
inline int InputDB::get(const std::string &key) const
{
  Require(check(key));
  return d_data_int.find(key)->second;
}
template <>
inline double InputDB::get(const std::string &key) const
{
  Require(check(key));
  return d_data_dbl.find(key)->second;
}
template <>
inline vec_int InputDB::get(const std::string &key) const
{
  Require(check(key));
  return d_data_vec_int.find(key)->second;
}
template <>
inline vec_dbl InputDB::get(const std::string &key) const
{
  Require(check(key));
  return d_data_vec_dbl.find(key)->second;
}
template <>
inline std::string InputDB::get(const std::string &key) const
{
  Require(check(key));
  return d_data_str.find(key)->second;
}

// Put

template <>
inline void InputDB::put(const std::string &key, const int value)
{
  d_data_int[key] = value;
  Ensure(check(key));
}
template <>
inline void InputDB::put(const std::string &key, const double value)
{
  d_data_dbl[key] = value;
  Ensure(check(key));
}
template <>
inline void InputDB::put(const std::string &key, const vec_int value)
{
  d_data_vec_int[key] = value;
  Ensure(check(key));
}
template <>
inline void InputDB::put(const std::string &key, const vec_dbl value)
{
  d_data_vec_dbl[key] = value;
  Ensure(check(key));
}
template <>
inline void InputDB::put(const std::string &key, const std::string value)
{
  d_data_str[key] = value;
  Ensure(check(key));
}

} // end namespace detran

#endif /* INPUTDB_I_HH_ */

//---------------------------------------------------------------------------//
//              end of InputDB.i.hh
//---------------------------------------------------------------------------//
