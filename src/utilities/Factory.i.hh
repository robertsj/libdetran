//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Factory.i.hh
 *  @brief Factory inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_utilities_FACTORY_I_HH_
#define detran_utilities_FACTORY_I_HH_

#include "DBC.hh"

namespace detran_utilities
{

//----------------------------------------------------------------------------//
template <class B>
Factory<B>& Factory<B>::Instance()
{
  static Factory instance;
  return instance;
}

//----------------------------------------------------------------------------//
template <class B>
bool Factory<B>::Register(const std::string &key, CreateFunction create)
{
  return Factory<B>::Instance().Register_private(key, create);
}

//----------------------------------------------------------------------------//
template <class B>
bool Factory<B>::IsRegistered(const std::string &key)
{
  return Factory<B>::Instance().d_callbacks.find(key)
         != Factory<B>::Instance().d_callbacks.end();
}

//----------------------------------------------------------------------------//
template <class B>
typename Factory<B>::CreateFunction
Factory<B>::GetCreateFunction(const std::string &key)
{
  typename CreateFunctionMap::const_iterator it = d_callbacks.find(key);
  if (it == d_callbacks.end())
  {
    THROW("Unknown key: " + key);
  }
  return it->second;
}

//----------------------------------------------------------------------------//
template <class B>
bool Factory<B>::Register_private(const std::string &key, CreateFunction create)
{
  return d_callbacks.insert(CreateFunctionMapPair(key, create)).second;
}

} // end namespace detran_utilities

#endif /* detran_utilities_FACTORY_I_HH_ */

//----------------------------------------------------------------------------//
//              end of Factory.i.hh
//----------------------------------------------------------------------------//

