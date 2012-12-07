//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   InputDB.hh
 *  @author Jeremy Roberts
 *  @brief  InputDB class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_utilities_INPUTDB_HH_
#define detran_utilities_INPUTDB_HH_

#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"
#include <string>
#include <map>
#ifdef DETRAN_ENABLE_BOOST
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/string.hpp>
#endif

namespace detran_utilities
{

//===========================================================================//
/**
 * \class InputDB
 * \brief Flexible storage for user input.
 *
 * User input for transport codes typically involves several integer and
 * scalar quantities along with meshing, materials data, and materials
 * placement.  For the former quantities, it can be a bonafide pain in the
 * arse to maintain an input structure throughout development.  To avoid that
 * issue, we use C++ maps of type int and double to store anything the user
 * needs to put there for use anywhere in the code.  This is inspired by the
 * way Denovo handles input.  For us, a descendent of Parser will read some
 * form of user input (be it XML or SILO or HDF5) and fill an InputDB along
 * with generating other things (like the data needed to make meshes and
 * materials libraries).  Of course, some parameters need to be there, and
 * that can be checked during a verification.
 *
 * (It's anticipated this will become a serment-wide approach)
 *
 */
//===========================================================================//

class InputDB
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  enum INPUT_TYPES
  {
    INT,
    DBL,
    STR,
    VEC_INT,
    VEC_DBL,
    SPINPUT
  };

  typedef SP<InputDB> SP_input;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor with optional name
  InputDB(std::string name = "InputDB");

  static SP_input Create(std::string name = "InputDB")
  {
    SP_input p;
    p = new InputDB(name);
    return p;
  }

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Return value of key
   *  @param    key     Name of the parameter.
   *  @param    value   Reference to which parameter value is assigned.
   *  @return           Check whether key is found.
   */
  template <class T>
  inline T get(const std::string &key) const
  {

  }

  /**
   *  @brief Check if the key is in a  map.
   *
   *  Note, if used consistently, this will prevent the
   *  same key being used for different value types.
   *  Maybe there's a better way to do this input
   *  structure.
   *
   *  @param    key     Name of the parameter.
   *  @return           True if exists; false otherwise.
   */
  inline bool check(const std::string &key) const;

  /**
   *  @brief Put a key and value in the database.
   *  @param    key         Name of the parameter.
   *  @param    value       Reference to which parameter value is assigned.
   *  @param    replace     Can we replace a current value?
   *  @return               Check whether key is found.
   */
  template <class T>
  inline void put(const std::string &key, const T value)
  {

  }

  /// Return a map
  template <class T>
  const std::map<std::string, T>& get_map();

  /// Number of entries of a certain type.
  int size(int type) const;

  /// Display all my contents.
  void display() const;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Name of the database
  std::string d_name;
  /// Integer parameters.
  std::map<std::string, int> d_data_int;
  /// Double parameters.
  std::map<std::string, double> d_data_dbl;
  /// Integer vector parameters.
  std::map<std::string, vec_int> d_data_vec_int;
  /// Double vector parameters.
  std::map<std::string, vec_dbl> d_data_vec_dbl;
  /// String parameters.
  std::map<std::string, std::string> d_data_str;
  /// InputDB parameters
  std::map<std::string, SP_input> d_data_db;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

#ifdef DETRAN_ENABLE_BOOST

  friend class boost::serialization::access;
  template<class Archive>

  void serialize(Archive & ar, const unsigned int version)
  {
    ar & d_data_int;
    ar & d_data_dbl;
    ar & d_data_vec_int;
    ar & d_data_vec_dbl;
    ar & d_data_str;
    ar & d_data_db;
  }

#endif

};

} // end namespace detran_utilities

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "InputDB.i.hh"

#endif /* detran_utilities_INPUTDB_HH_ */

//---------------------------------------------------------------------------//
//              end of InputDB.hh
//---------------------------------------------------------------------------//
