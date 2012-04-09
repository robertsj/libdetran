//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GenDB.hh
 * \author Jeremy Roberts
 * \date   Mar 19, 2012
 * \brief  
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
// $Rev::                                               $:Rev of last commit
// $Author:: j.alyn.roberts@gmail.com                   $:Author of last commit
// $Date::                                              $:Date of last commit
//---------------------------------------------------------------------------//


#ifndef GENDB_HH_
#define GENDB_HH_

#include <string>
#include <map>

#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

class GenDB : public Object
{

public:

  /*
   *  \brief Constructor.
   */
  GenDB();

  /// \name Accessors
  //\{

  /*
   *  \brief Return value of key by reference.
   *  \param    key     Name of the parameter.
   *  \param    value   Reference to which parameter value is assigned.
   *  \return           Check whether key is found.
   */
//  template <class T>
//  bool get(const std::string &key, T &value) const;

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

  //\}

  /*
   *  \brief Validate the database.
   *
   *  This is probably the only place where upkeep needs to occur.  There
   *  may be an intelligent way to do it, something like a "grammar".
   *
   *  \return Whether or not verification was successful.
   */
  bool is_valid() const {};

private:

  /*
   *  \brief Empty map wrapper
   */
  struct BaseItem
  {
  };

  /*
   *  \brief Template map wrapper
   */
  template <class T>
  struct TemplateItem : public BaseItem
  {
    explicit TemplateItem(T item){d_item = item;}
    T d_item;
  };

  typedef SP<BaseItem> SP_baseitem;
  typedef std::map<std::string, SP_baseitem> map_type;
  typedef std::pair<std::string, SP_baseitem> map_pair;

  /// \name Data
  //\{

  /// Vector of base map pointers.
  map_type d_map;

  //\}

  /// \name Implementation
  //\{

  //\}

};

} // end namespace detran

#endif /* GENDB_HH_ */

//---------------------------------------------------------------------------//
//              end of GenDB.hh
//---------------------------------------------------------------------------//
