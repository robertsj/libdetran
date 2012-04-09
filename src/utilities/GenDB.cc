//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GenDB.cc
 * \author Jeremy Roberts
 * \date   Mar 19, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "GenDB.hh"

// System
#include <iostream>

namespace detran
{

GenDB::GenDB()
{
}

template <class T>
void GenDB::put(const std::string &key, const T value)
{
  // Erase the value associated with the key if it exists
  map_type::iterator it;
  it = d_map.find(key);
  if (it != d_map.end())
    d_map.erase(it);

  // Make the new value and add it.
  SP_baseitem sp_val;
  sp_val = new TemplateItem<T> (value);
  d_map[key] = sp_val;
}

template <>
void GenDB::put(const std::string &key, const int value)
{
  // Make the new value and add it.
  SP_baseitem sp_val;
  sp_val = new TemplateItem<int> (value);
  d_map[key] = sp_val;
}

// I'm not sure if map.end()->second has a usable value...
template <class T>
T GenDB::get(const std::string &key) const
{
  typedef SP<TemplateItem<T> > SP_titem;
  SP_titem sp_val = (SP_titem) d_map.find(key)->second;
  return *sp_val;
}

//template <>
//int GetDB::get(const std::string &key) const
//{
//  typedef SP<TemplateItem<int> > SP_titem;
//  SP_titem sp_val = (SP_titem) d_map.find(key)->second;
//  return *sp_val;
//}

} // end namespace detran


//---------------------------------------------------------------------------//
//              end of GenDB.cc
//---------------------------------------------------------------------------//
