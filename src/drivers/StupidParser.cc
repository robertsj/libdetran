//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   StupidParser.cc
 *  @brief  StupidParser member definitions
 *  @note   Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#include "geometry/Mesh1D.hh"
#include "geometry/Mesh2D.hh"
#include "geometry/Mesh3D.hh"
#include "drivers/StupidParser.hh"
#include <sstream>
#include <iostream>
#include <ios>

namespace detran
{

//---------------------------------------------------------------------------//
StupidParser::StupidParser(int argc,  char **argv)
  : d_is_hdf5(false)
{
  Insist(argc > 1, "Not enough command line arguments!");
  d_filename = std::string(argv[1]);
  open_file();
}

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
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

//---------------------------------------------------------------------------//
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
detran_geometry::Mesh::SP_mesh StupidParser::parse_mesh_text()
{
  Insist(d_input, "Input must be parsed before mesh!");
  using namespace detran_geometry;

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
      Insist(d_input->check("mesh_ycme"),
             "mesh_xcme given, but not ycme for a 2D problem!");
      yme = d_input->get<vec_dbl>("mesh_ycme");
      Insist(d_input->check("mesh_yfm"),
             "mesh_ycme given, but not yfm!");
      yfm = d_input->get<vec_int>("mesh_yfm");
      Insist(xme.size() == yfm.size() + 1,
             "size(ycme) must be size(yfm)+1");
      mesh_map_size *= yfm.size();
    }
    else
    {
      Insist(d_input->check("mesh_xfme"),
             "mesh_xfme given, but not yfme for a 2D problem!");
      yme = d_input->get<vec_dbl>("mesh_yfme");
      mesh_map_size *= (yme.size()-1);
    }
  }
  if (dimension > 2)
  {
    if (coarse)
    {
      Insist(d_input->check("mesh_zcme"),
             "mesh_xcme given, but not zcme for a 3D problem!");
      zme = d_input->get<vec_dbl>("mesh_zcme");
      Insist(d_input->check("mesh_zfm"),
             "mesh_zcme given, but not zfm!");
      zfm = d_input->get<vec_int>("mesh_zfm");
      Insist(zme.size() == zfm.size() + 1,
             "size(zcme) must be size(zfm)+1");
      mesh_map_size *= zfm.size();
    }
    else
    {
      Insist(d_input->check("mesh_zfme"),
             "mesh_xfme given, but not zfme for a 3D problem!");
      zme = d_input->get<vec_dbl>("mesh_zfme");
      mesh_map_size *= (zme.size()-1);
    }
  }
  Insist(d_input->check("mesh_mat_map"), "mesh_mat_map must be defined!");
  mesh_map = d_input->get<vec_int>("mesh_mat_map");
  Insist(mesh_map.size() == mesh_map_size, "mesh_mat_map has wrong size!");

  if (dimension == 1)
  {
    if (coarse)  d_mesh = new Mesh1D(xfm, xme, mesh_map);
    if (!coarse) d_mesh = new Mesh1D(xme, mesh_map);
  }
  else if (dimension == 2 && coarse)
  {
    if (coarse)  d_mesh = new Mesh2D(xfm, yfm, xme, yme, mesh_map);
    if (!coarse) d_mesh = new Mesh2D(xme, yme, mesh_map);
  }
  else if (dimension == 3 && coarse)
  {
    if (coarse)  d_mesh = new Mesh3D(yfm, zfm, zfm, xme, yme, zme, mesh_map);
    if (!coarse) d_mesh = new Mesh3D(xme, yme, zme, mesh_map);
  }
  else
  {
    THROW("Something wierd is going on here...");
  }

  return d_mesh;
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of StupidParser.cc
//---------------------------------------------------------------------------//



