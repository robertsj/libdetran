//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SimpleDB.hh
 * \author Jeremy Roberts
 * \brief  SimpleDB class definition.
 */
//---------------------------------------------------------------------------//

#ifndef SIMPLEDB_HH_
#define SIMPLEDB_HH_

#include <string>
#include <map>

#include "DBC.hh"
#include "Definitions.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class SimpleDB
 * \brief Relatively flexible storage for user input.
 *
 * I think as implemented, this needs to be a singleton class!
 *
 */

class SimpleDB : public Object
{

public:

  /*
   *  \brief Constructor.
   */
  SimpleDB();

  /// \name Accessors
  //\{

  /*
   *  \brief Return value of key
   *  \param    key     Name of the parameter.
   *  \param    value   Reference to which parameter value is assigned.
   *  \return           Check whether key is found.
   */
  template <class T>
  void get(std::string key, T value) const;

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
  struct BaseMap
  {
  };

  /*
   *  \brief Template map wrapper
   */
  template <class T>
  struct TemplateMap : public BaseMap
  {
    std::map<std::string, T> d_data;
  };
  
  /// Smart pointer to empty map wrapper.
 // typedef typename SP<BaseMap> SP_basemap;

  template<class T>
  struct type_index
  {
    static int val;
    static int exists()
    {
      if (val == -1)
      {
        val = ++SimpleDB::d_count;
        return 0;
      }
      return 1;
    }
  };


  /// \name Data
  //\{

  /// Count of mapped types.
  static int d_count;

  /// Vector of base map pointers.
  std::vector<BaseMap*>  d_maps;

  //\}

  /// \name Implementation
  //\{

  //\}

};

// Initialize type counter.
int SimpleDB::d_count = -1;

// Initialize type index.
template <typename T>
int SimpleDB::type_index<T>::val = -1;


} // end namespace detran

#endif /* SIMPLEDP_HH_ */

//---------------------------------------------------------------------------//
//              end of SimpleDB.hh
//---------------------------------------------------------------------------//
