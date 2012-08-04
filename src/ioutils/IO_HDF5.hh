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

// Configuration
#include "detran_config.h"

// Detran
#include "Material.hh"
#include "Mesh.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"
#include "SP.hh"

// System
#include <string>

#ifdef DETRAN_ENABLE_HDF5
#include "hdf5.h"
#include "IO_HDF5_Traits.hh"
#endif

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

  typedef detran::SP<IO_HDF5>             SP_io_hdf5;
  typedef detran::InputDB::SP_input       SP_input;
  typedef detran::Material::SP_material   SP_material;
  typedef detran::Mesh::SP_mesh           SP_mesh;
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
   *  \brief Write the material database into an HDF5 file.
   *  \param mat    Material database to be written
   */
  void write(SP_material mat);

  /*!
   *  \brief Write the material database into an HDF5 file.
   *  \param mat    Material database to be written
   */
  void write(SP_mesh mesh);

  /// Close an HDF5 file if open.
  void close();

  // Get an input database from file.
  SP_input read_input();

  /// Get a material database from file.
  SP_material read_material();

  /// Get a mesh from file.
  SP_mesh read_mesh();

  bool is_valid() const
  {
    return true;
  }

private:

#ifdef DETRAN_ENABLE_HDF5

  /// \name Private Data
  /// \{

  /// HDF5 file id
  hid_t d_file_id;

  /// HDF5 filename
  std::string d_filename;

  /// HDF5 is open
  bool d_open;

  SP_input d_input;
  SP_material d_material;
  SP_mesh d_mesh;

  /// \}

  /// \name Implementation
  /// \{

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

  /// Write a map to HDF5.
  template <class T>
  bool write_map(hid_t group,
                 const char* name,
                 const std::map<std::string, T> &map);

  /// Read a map from HDF5.
  template <class T>
  bool read_map(hid_t group,
                const char* name,
                std::map<std::string, T> &map);

  /// Write into a vector (int or double) to HDF5.
  template <class T>
  bool write_vec(hid_t group, const char* name, const std::vector<T> &source);

  /// Read into a vector (int or double)
  template <class T>
  bool read_vec(hid_t group, const char* name, std::vector<T> &target);

  /// Write a scalar (int or double) attribute
  template <class T>
  bool write_scalar_attribute(hid_t group, const char* name, const T &value);

  /// Read a scalar (int or double) attribute
  template <class T>
  bool read_scalar_attribute(hid_t group, const char* name, T &value);

  /// Check if a location (group, dataset, etc.) "name" exists
  bool exists(hid_t location, const char* name)
  {
    htri_t flag = H5Lexists(location, name, H5P_DEFAULT);
    if (flag <= 0) return false;
    return true;
  }

  /// \}

#endif

};

} // end namespace detran_ioutils

//---------------------------------------------------------------------------//
// TEMPLATE DEFINITIONS
//---------------------------------------------------------------------------//

#ifdef DETRAN_ENABLE_HDF5
#include "IO_HDF5.t.hh"
#endif

#endif // IO_HDF5_HH_

//---------------------------------------------------------------------------//
//              end of file IO_HDF5.hh
//---------------------------------------------------------------------------//
