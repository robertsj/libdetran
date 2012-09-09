//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   StupidParser.cc
 * \author robertsj
 * \date   Apr 11, 2012
 * \brief  StupidParser member definitions.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "Mesh1D.hh"
#include "Mesh2D.hh"
#include "Mesh3D.hh"

// Utilities
#include "StupidParser.hh"
#include "Definitions.hh"

// System
#include <sstream>
#include <iostream>
#include <ios>

namespace detran
{

StupidParser::StupidParser(int argc,  char **argv)
  : d_is_hdf5(false)
{
  Insist(argc > 1, "Not enough command line arguments!");

  std::string filename(argv[1]);
  std::string ext = get_file_extension(filename);

  if (ext == "h5")
  {
#ifdef DETRAN_ENABLE_HDF5
    // Open the HDF5 file.
    d_hdf5 = new detran_ioutils::IO_HDF5(filename);
    d_is_hdf5 = true;
#else
    THROW("User specified HDF5 input but HDF5 is not enabled!");
#endif
  }
  else
  {
    // Open the file.
    d_file.open(argv[1]);
    if (!d_file) THROW("The file could not be opened!");
  }

}

StupidParser::SP_input StupidParser::parse_input()
{
  SP_input input;
  // Get the input
  if (!d_is_hdf5)
    input = parse_input_text();
  else
  {
#ifdef DETRAN_ENABLE_HDF5
    input = d_hdf5->read_input();
#endif
  }
  // Postconditions
  Insist(input, "Error parsing input.")
  return input;
}

StupidParser::SP_material StupidParser::parse_material()
{
  // Get the material
  SP_material mat;
  if (!d_is_hdf5)
    mat = parse_material_text();
  else
  {
#ifdef DETRAN_ENABLE_HDF5
    mat = d_hdf5->read_material();
#endif
  }
  // Postconditions
  Insist(mat, "Error parsing material.")
  return mat;
}

StupidParser::SP_mesh StupidParser::parse_mesh()
{
  // Get the mesh
  SP_mesh mesh;
  if (!d_is_hdf5)
    mesh = parse_mesh_text();
  else
  {
#ifdef DETRAN_ENABLE_HDF5
    mesh = d_hdf5->read_mesh();
#endif
  }
  // Postconditions
  Insist(mesh, "Error parsing mesh.")
  return mesh;
}

//---------------------------------------------------------------------------//
// IMPLEMENTATION
//---------------------------------------------------------------------------//

/// Parse input. \todo maybe template some of this?
StupidParser::SP_input StupidParser::parse_input_text()
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

/// Parse material.
detran_material::Material::SP_material
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
          downscatter = read_scalar<int>(iss);
        }
      }
    }
  }
  Require(number_materials > 0);

  d_material = new detran_material::
    Material(number_materials, number_groups, downscatter);

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
          Insist(value.size() == number_groups, "diff_coef coefficient wrong length.");
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

detran_geometry::Mesh::SP_mesh
StupidParser::parse_mesh_text()
{
  Insist(d_input, "Input must be parsed before mesh!");

  bool coarse = false;
  Insist(d_input->check("dimension"), "The dimension must be specified!");
  int dimension = d_input->get<int>("dimension");

  vec_dbl xme;
  vec_dbl yme;
  vec_dbl zme;
  vec_int xfm;
  vec_int yfm;
  vec_int zfm;
  vec_int mesh_map;
  int mesh_map_size = 1;

  if (d_input->check("mesh_xcme"))
  {
    xme = d_input->get<vec_dbl>("mesh_xcme");
    coarse = true;
    Insist(d_input->check("mesh_xfm"), "mesh_xcme given, but not xfm!");
    xfm = d_input->get<vec_int>("mesh_xfm");
    Insist(xme.size() == xfm.size() + 1, "size(xcme) must be size(xfm)+1");
    mesh_map_size *= xfm.size();
  }
  else
  {
    Insist(d_input->check("mesh_xfme"), "either xcme or xfme must be defined!");
    xme = d_input->get<vec_dbl>("mesh_xfme");
    mesh_map_size *= (xme.size()-1);
  }

  if (dimension > 1)
  {
    if (coarse)
    {
      Insist(d_input->check("mesh_ycme"), "mesh_xcme given, but not ycme for a 2D problem!");
      yme = d_input->get<vec_dbl>("mesh_ycme");
      Insist(d_input->check("mesh_yfm"), "mesh_ycme given, but not yfm!");
      yfm = d_input->get<vec_int>("mesh_yfm");
      Insist(xme.size() == yfm.size() + 1, "size(ycme) must be size(yfm)+1");
      mesh_map_size *= yfm.size();
    }
    else
    {
      Insist(d_input->check("mesh_xfme"), "mesh_xfme given, but not yfme for a 2D problem!");
      yme = d_input->get<vec_dbl>("mesh_yfme");
      mesh_map_size *= (yme.size()-1);
    }
  }
  if (dimension > 2)
  {
    if (coarse)
    {
      Insist(d_input->check("mesh_zcme"), "mesh_xcme given, but not zcme for a 3D problem!");
      zme = d_input->get<vec_dbl>("mesh_zcme");
      Insist(d_input->check("mesh_zfm"), "mesh_zcme given, but not zfm!");
      zfm = d_input->get<vec_int>("mesh_zfm");
      Insist(zme.size() == zfm.size() + 1, "size(zcme) must be size(zfm)+1");
      mesh_map_size *= zfm.size();
    }
    else
    {
      Insist(d_input->check("mesh_zfme"), "mesh_xfme given, but not zfme for a 3D problem!");
      zme = d_input->get<vec_dbl>("mesh_zfme");
      mesh_map_size *= (zme.size()-1);
    }
  }
  Insist(d_input->check("mesh_mat_map"), "mesh_mat_map must be defined!");
  mesh_map = d_input->get<vec_int>("mesh_mat_map");
  Insist(mesh_map.size() == mesh_map_size, "mesh_mat_map has wrong size!");

  if (dimension == 1)
  {
    if (coarse)  d_mesh = new detran_geometry::Mesh1D(xfm, xme, mesh_map);
    if (!coarse) d_mesh = new detran_geometry::Mesh1D(xme, mesh_map);
  }
  else if (dimension == 2 and coarse)
  {
    if (coarse)  d_mesh = new detran_geometry::Mesh2D(xfm, yfm, xme, yme, mesh_map);
    if (!coarse) d_mesh = new detran_geometry::Mesh2D(xme, yme, mesh_map);
  }
  else if (dimension == 3 and coarse)
  {
    if (coarse)  d_mesh = new detran_geometry::Mesh3D(yfm, zfm, zfm, xme, yme, zme, mesh_map);
    if (!coarse) d_mesh = new detran_geometry::Mesh3D(xme, yme, zme, mesh_map);
  }
  else
  {
    THROW("Something wierd is going on here...");
  }

  return d_mesh;
}

//
void StupidParser::remove_white_space(std::string &s)
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

//
template <class T>
std::vector<T>
StupidParser::read_vector(std::istringstream &iss)
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
//
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

//---------------------------------------------------------------------------//
//              end of StupidParser.cc
//---------------------------------------------------------------------------//



