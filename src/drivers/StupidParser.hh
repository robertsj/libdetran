//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   StupidParser.hh
 *  @brief  StupidParser class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_STUPIDPARSER_HH_
#define detran_STUPIDPARSER_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "external_source/ExternalSource.hh"
#ifdef DETRAN_ENABLE_HDF5
#include "ioutils/IO_HDF5.hh"
#endif
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class StupidParser
 * \brief Parse a really simple input file.
 *
 * Because I'm exceedingly lazy and want to use Python for all my "input"
 * now that I've tried it, I really want to spend little time on a text-based
 * input format.  Hence, I'm doing the absolutely easiest (IMHO) thing
 * possible.  For those who need more flexibility, write a wrapper or
 * better, contribute a more reasonable parser.
 * The input file will be a series of lines similar to the
 * following, which is a simple working example
 *
 * \verbatim

   # General input (the "#" is ignored)
   int      number_groups           2
   int      dimension               2
   dbl      inner_max_tolerance     1e-4
   str      problem_type            fixed

   # Mesh (blank lines also ignored)
   mesh mesh_xcme 0.0 1.0 2.0
   mesh mesh_ycme 0.0 2.0 4.0
   mesh mesh_xfm  10
   mesh mesh_yfm  20
   # all arrays are read in 1d; logically these ordered by x->y->z
   mesh mesh_mat  0 0 1 0

   # Material
   material number_materials 2
   # moderator
   material sigma_t    0   0.1890 1.4633
   material sigma_s    0 0 0.1507 0.0000
   material sigma_s    0 1 0.0380 1.4536
   # fuel
   material sigma_t    1   0.2263 1.0119
   material sigma_s    1 0 0.2006 0.0000
   material sigma_s    1 1 0.0161 0.9355
   material nu_sigma_f 1   0.0067 0.1241
   material chi        1   1.0000 0.0000

   # Fixed source (only isotropic or constant)
   source type isotropic
   source number_spectra 2
   source spectrum 0 1.0 0.0
   source spectrum 1 0.0 0.0
   source map      0 0 1 0

 * \endverbatim
 *
 * \sa InputDB
 */
//---------------------------------------------------------------------------//

class StupidParser
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_material::Material::SP_material  SP_material;
#ifdef DETRAN_ENABLE_HDF5
  typedef detran_ioutils::IO_HDF5::SP_io_hdf5     SP_io_hdf5;
#endif
  typedef detran_utilities::vec_int               vec_int;
  typedef detran_utilities::vec_dbl               vec_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  @param    argc    Number of arguments
   *  @param    argv    Arguments
   */
  StupidParser(int argc, char **argv);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Parse the input.
  SP_input parse_input();

  /// Parse mesh.
  SP_mesh parse_mesh();

  /// Parse material.
  SP_material parse_material();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_input d_input;
  SP_mesh d_mesh;
  SP_material d_material;

  /// File handle
  std::ifstream d_file;

  /// HDF5 file handler.
#ifdef DETRAN_ENABLE_HDF5
  SP_io_hdf5 d_hdf5;
#endif

  /// HDF5 enabled flag
  bool d_is_hdf5;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Get input from text file
  SP_input parse_input_text();

  /// Get material from text file
  SP_material parse_material_text();

  /// Get mesh from text file
  SP_mesh parse_mesh_text();

  /// Remove double spaces
  void remove_white_space(std::string &s);

  template <class T>
  T read_scalar(std::istringstream &iss);

  template <class T>
  std::vector<T> read_vector(std::istringstream &iss);

  /// Get the file extension
  std::string get_file_extension(const std::string& filename)
  {
      if(filename.find_last_of(".") != std::string::npos)
        return filename.substr(filename.find_last_of(".") + 1);
      return "";
  }

};

} // end namespace detran

#endif /* detran_STUPIDPARSER_HH_ */

//---------------------------------------------------------------------------//
//              end of StupidParser.hh
//---------------------------------------------------------------------------//
