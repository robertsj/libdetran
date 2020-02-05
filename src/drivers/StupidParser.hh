//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file   StupidParser.hh
 *  @brief  StupidParser class definition
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

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

//----------------------------------------------------------------------------//
/**
 *  @class StupidParser
 *  @brief Parse a really simple input file.
 *
 *  Because I'm exceedingly lazy and want to use Python for all my "input"
 *  now that I've tried it, I really want to spend little time on a text-based
 *  input format.  Hence, I'm doing the absolutely easiest (IMHO) thing
 *  possible.  For those who need more flexibility, write a wrapper or
 *  better, contribute a more reasonable parser.
 *
 *  The input file will be a series of lines similar to the
 *  following, which is a simple working example.  Note, all lines containing
 *  data must be flush left.
 *
 *  @verbatim

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

 *  @endverbatim
 *
 *  @sa InputDB
 */
//----------------------------------------------------------------------------//

class StupidParser
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input     SP_input;
  typedef detran_geometry::Mesh::SP_mesh          SP_mesh;
  typedef detran_material::Material::SP_material  SP_material;
#ifdef DETRAN_ENABLE_HDF5
  typedef detran_ioutils::IO_HDF5::SP_io_hdf5     SP_io_hdf5;
#endif
  typedef detran_utilities::vec_int               vec_int;
  typedef detran_utilities::vec_dbl               vec_dbl;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param    argc    Number of arguments
   *  @param    argv    Arguments
   */
  StupidParser(int argc, char **argv);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Parse the input.
  SP_input parse_input();

  /// Parse mesh.
  SP_mesh parse_mesh();

  /// Parse material.
  SP_material parse_material();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  SP_input d_input;
  SP_mesh d_mesh;
  SP_material d_material;
  std::string d_filename;
  /// File handle
  std::ifstream d_file;

  /// HDF5 file handler.
#ifdef DETRAN_ENABLE_HDF5
  SP_io_hdf5 d_hdf5;
#endif

  /// HDF5 enabled flag
  bool d_is_hdf5;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

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

  void open_file()
  {
	std::cout << "Trying to open user file: " << d_filename << std::endl;

	std::string ext = "";
	if (d_filename.length() > 3)
	  ext = d_filename.substr(d_filename.length()-2, 2);
    // if(filename.find_last_of(".") < filename.size())
    //   ext = filename.substr(filename.find_last_of(".") + 1);

    if (ext == "h5")
    {
#ifdef DETRAN_ENABLE_HDF5
      // Open the HDF5 file.
      d_hdf5 = new detran_ioutils::IO_HDF5(d_filename);
      d_is_hdf5 = true;
#else
      THROW("User specified HDF5 input but HDF5 is not enabled!");
#endif
    }
    else
    {
      // Open the file.
      d_file.open(d_filename.c_str());
      if (!d_file) THROW("The file could not be opened!");
    }
	std::cout << "The file was opened successfully." << std::endl;
  }

};

//----------------------------------------------------------------------------//
inline StupidParser::SP_input 
StupidParser::parse_input_text()
{
  using std::cout;
  using std::endl;
  using std::string;
  using std::getline;
  using std::istringstream;

  Insist(d_file.is_open(), "The input file must be open to parse!");

  SP_input p;
  p = new detran_utilities::InputDB();

  string line;
  while (getline(d_file, line))
  {
    remove_white_space(line);
#ifdef DETRAN_ENABLE_DEBUG
    // dump reformatted input back out
    cout << line << endl;
#endif
    istringstream iss(line);

    // Skip comments
    if (iss.peek() == '#') continue;

    string word;
    string key;
    if (getline(iss, word, ' '))
    {
      if (word == "int")
      {
        key = read_scalar<string>(iss);
        Require(key.size() > 0);
        p->put<int>(key, read_scalar<int>(iss));
      }
      else if (word == "dbl")
      {
        key = read_scalar<string>(iss);
        Require(key.size() > 0);
        p->put<double>(key, read_scalar<double>(iss));
      }
      else if (word == "str")
      {
        key = read_scalar<string>(iss);
        Require(key.size() > 0);
        p->put<string>(key, read_scalar<string>(iss));
      }
      else if (word == "vec_int")
      {
        key = read_scalar<string>(iss);
        Require(key.size() > 0);
        p->put<vec_int>(key, read_vector<int>(iss));
      }
      else if (word == "vec_dbl")
      {
        key = read_scalar<string>(iss);
        Require(key.size() > 0);
        p->put<vec_dbl>(key, read_vector<double>(iss));
      }
    }
  }
  d_input = p;
  return p;
}

//----------------------------------------------------------------------------//
inline detran_material::Material::SP_material 
StupidParser::parse_material_text()
{
  using std::cout;
  using std::endl;
  using std::string;
  using std::getline;
  using std::istringstream;

  Insist(d_file.is_open(), "The input file must be open to parse!");
  Insist(d_input, "Input must be parsed before material!");

  d_file.clear();
  d_file.seekg(0, std::ios::beg);

  string line;

  Require(d_input->check("number_groups"));
  int number_groups = d_input->get<int>("number_groups");
  int number_materials = 0;
  bool downscatter = false;

  // Go through and get number of materials.
  while (getline(d_file, line))
  {
    remove_white_space(line);
    istringstream iss(line);
    if (iss.peek() == '#') continue;
    string val;
    if (getline(iss, val, ' '))
    {
      if (val == "material")
      {

        val = read_scalar<string>(iss);
        if (val == "number_materials")
        {
          number_materials = read_scalar<int>(iss);
        }
        else if (val == "downscatter")
        {
          downscatter = 0 != read_scalar<int>(iss);
        }
      }
    }
  }
  Require(number_materials > 0);

  d_material = new detran_material::Material(number_materials, number_groups);
  d_material->set_downscatter(downscatter);

  // Fill the rest.
  d_file.clear();
  d_file.seekg(0, std::ios::beg);
  while (getline(d_file, line))
  {
    remove_white_space(line);
    istringstream iss(line);
    if (iss.peek() == '#') continue;
    string val;
    if (getline(iss, val, ' '))
    {
      if (val == "material")
      {
        val = read_scalar<string>(iss);
        if (val == "sigma_t")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "sigma_t wrong length.");
          d_material->set_sigma_t(mat, value);
        }
        if (val == "sigma_a")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "sigma_a wrong length.");
          d_material->set_sigma_a(mat, value);
        }
        if (val == "sigma_f")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "sigma_f wrong length.");
          d_material->set_sigma_f(mat, value);
        }
        if (val == "nu")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "nu wrong length.");
          d_material->set_nu(mat, value);
        }
        if (val == "chi")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "chi wrong length.");
          d_material->set_chi(mat, value);
        }
        if (val == "diff_coef")
        {
          int mat = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups,
                 "diff_coef coefficient wrong length.");
          d_material->set_diff_coef(mat, value);
        }
        if (val == "sigma_s")
        {
          int mat = read_scalar<int>(iss);
          int g   = read_scalar<int>(iss);
          vec_dbl value = read_vector<double>(iss);
          Insist(value.size() == number_groups, "sigma_s wrong length.");
          d_material->set_sigma_s(mat, g, value);
        }
      }
    }
  } // end material loop
  d_material->finalize();
  return d_material;
}

//----------------------------------------------------------------------------//
inline void StupidParser::remove_white_space(std::string &s)
{
  int i = s.find("  ");
  if (i > -1)
  {
    s.replace(i, 2, " ");
    remove_white_space(s);
  }
  else if (i == -1)
  {
    return;
  }
}

//----------------------------------------------------------------------------//
template <class T>
std::vector<T> StupidParser::read_vector(std::istringstream &iss)
{
  std::vector<T> vec;
  std::string s;
  while (std::getline(iss, s, ' '))
  {
    double v;
    std::istringstream tmp;
    tmp.str(s);
    tmp >> v;
    vec.push_back(v);
  }
  return vec;
}

//----------------------------------------------------------------------------//
template <class T>
T StupidParser::read_scalar(std::istringstream &iss)
{
  T val;
  std::string s;
  getline(iss, s, ' ');
  std::istringstream tmp;
  tmp.str(s);
  tmp >> val;
  return val;
}

} // end namespace detran

#endif /* detran_STUPIDPARSER_HH_ */

//----------------------------------------------------------------------------//
//              end of StupidParser.hh
//----------------------------------------------------------------------------//
