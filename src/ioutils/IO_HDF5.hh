//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   IO_HDF5.hh
 * \brief  IO_HDF5
 * \author Jeremy Roberts
 * \date   Jul 29, 2012
 */
//---------------------------------------------------------------------------//

#ifndef IO_HDF5_HH_
#define IO_HDF5_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "IO_HDF5_Traits.hh"

// System
#include <string>
#include "hdf5.h"

namespace detran_ioutils
{

/*!
 *  \class IO_HDF5
 *  \brief Maps an HDF5 database to an InputDB and vice versa.
 *
 *  A problem can be specified completely using the input,
 *  a material definition, and the mesh definition, all of
 *  which can actually live in the InputDB (as done in
 *  \ref StupidParser).  By mapping an HDF5 file to an
 *  input object, we can eliminate the fickle text processing,
 *  using Python or another interface to construct the HDF5
 *  file.  Then, the executable detran can be used, which is
 *  really handy for profiling and other analyses.
 *
 *  For now, we'll use a single group, /input.  All
 *  entries will go in that group.  Later, it might be useful
 *  to add ones for the mesh and material specification.
 */
class IO_HDF5: public detran::Object
{

public:

  typedef detran::InputDB::SP_input       SP_input;
  typedef detran::Material::SP_material   SP_material;
  typedef detran::vec_int                 vec_int;
  typedef detran::vec_dbl                 vec_dbl;

  /*!
   *  \brief Constructor
   *  \param filename HDF5 filename
   */
  IO_HDF5(std::string filename);

  /// Open HDF5 file for writing.  This replaces old content.
  void open();

  /*!
   *  \brief Write the input database into an HDF5 file.
   *  \param input    Input database to be written
   */
  void write(SP_input input);

  /*!
   *  \brief Write the input database into an HDF5 file.
   *  \param input    Input database to be written
   */
  void write(SP_material material);

  /// Close an HDF5 file if open.
  void close();

  /*!
   *  \brief Fill an input database from an HDF5 file.
   *  \param input    Input database to be filled
   *  \param filename HDF5 filename
   */
  void read(SP_input input);

  bool is_valid() const
  {
    return true;
  }

private:


  /// \name Private Data
  /// \{

  /// HDF5 file id
  hid_t   d_file_id;

  /// HDF5 filename
  std::string d_filename;

  /// Variable string type
  hid_t d_string_type;

  /// HDF5 is open
  bool d_open;

  /// \}

  /// \name Implementation
  /// \{

  /*!
   *  \brief Fill a temporary compound type container and write to file
   *  \param input  User input database
   *  \param data   Pointer to compound type array
   */
  template <class T>
  bool write_data(SP_input input,
                  hid_t group,
                  std::string name,
                  compound_type<T> *data);

  /*!
   *  \brief Fill a temporary compound type container and write to file
   *  \param input  User input database
   *  \param data   Pointer to compound type array
   */
  template <class T>
  bool read_data(SP_input input, hid_t group, std::string name);

  /// Set the data type for storage in memory
  template <class T>
  hid_t set_memtype();

  /// Set the data type for storage on disk
  template <class T>
  hid_t set_filetype();

  /// \}

};

} // end namespace detran_ioutils

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#include "IO_HDF5.t.hh"

#endif // IO_HDF5_HH_

//---------------------------------------------------------------------------//
//              end of file IO_HDF5.hh
//---------------------------------------------------------------------------//
