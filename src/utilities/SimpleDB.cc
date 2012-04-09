//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   InputDB.c
 * \author Jeremy Roberts
 * \brief  SimpleDB class member definition.
 */
//---------------------------------------------------------------------------//


#include "SimpleDB.hh"

#include <iostream>

namespace detran
{

SimpleDB::SimpleDB()
{
}


template <class T>
void SimpleDB::get(std::string key, T value) const
{
  if (type_index<T>::exists() == 0)
    std::cout << "This type doesn't have a map." << std::endl;
  // Return the map
  TemplateMap<T> *Tmap = (TemplateMap<T> *) d_maps[type_index<T>::val];
  value = Tmap->d_data.find(key)->second;
}


template <class T>
void SimpleDB::put(const std::string &key, const T value)
{
  if (type_index<T>::exists() == 0)
  { 
    //std::cout << " ADDING: " << type << std::endl; 
    BaseMap *m;
    m = new TemplateMap<T>;
    //m->me = type_index<T>::val;
    d_maps.push_back(m);
  }
  // Return the map
  std::cout << " Tval = " << type_index<T>::val << std::endl;
  std::cout << " vec size is " << d_maps.size() << std::endl;
  TemplateMap<T> *Tmap = (TemplateMap<T> *) d_maps[type_index<T>::val];
  // BaseMap *m = _maps[Tidx<T>::val];
  Tmap->d_data[key] = value;
}


} // end namespace detran

//---------------------------------------------------------------------------//
//              end of SimpleDB.cc
//---------------------------------------------------------------------------//
