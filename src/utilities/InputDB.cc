//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputDB.c
 * \author Jeremy Roberts
 * \brief  InputDB class definition.
 */
//---------------------------------------------------------------------------//


#include "InputDB.hh"

#include <iostream>

namespace detran_utils
{

InputDB::InputDB()
{
}


template <class T>
T InputDB::get(const std::string &key) const
{
}
template <>
int InputDB::get(const std::string &key) const
{
  return d_data_int.find(key)->second;
}
template <>
double InputDB::get(const std::string &key) const
{
  return d_data_dbl.find(key)->second;
}
template <>
vec_int InputDB::get(const std::string &key) const
{
  return d_data_vec_int.find(key)->second;
}
template <>
vec_dbl InputDB::get(const std::string &key) const
{
  return d_data_vec_dbl.find(key)->second;
}
template <>
std::string InputDB::get(const std::string &key) const
{
  return d_data_str.find(key)->second;
}


template <class T>
void InputDB::put(const std::string &key, const T value)
{
}
template <>
void InputDB::put(const std::string &key, const int value)
{
  d_data_int[key] = value;
}
template <>
void InputDB::put(const std::string &key, const double value)
{
  d_data_dbl[key] = value;
}
template <>
void InputDB::put(const std::string &key, const vec_int value)
{
  d_data_vec_int[key] = value;
}
template <>
void InputDB::put(const std::string &key, const vec_dbl value)
{
  d_data_vec_dbl[key] = value;
}
template <>
void InputDB::put(const std::string &key, const std::string value)
{
  d_data_str[key] = value;
}

} // end namespace detran_utils

//---------------------------------------------------------------------------//
//              end of InputDB.cc
//---------------------------------------------------------------------------//
