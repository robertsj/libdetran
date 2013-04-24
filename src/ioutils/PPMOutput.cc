//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PPMOutput.cc
 *  @author robertsj
 *  @date   Mar 13, 2013
 *  @brief  PPMOutput class definition.
 */
//---------------------------------------------------------------------------//

#include "ioutils/PPMOutput.hh"
#include "utilities/MathUtilities.hh"
#include <fstream>
#include <cstdio>
#include <sstream>

namespace detran_ioutils
{

//---------------------------------------------------------------------------//
PPMOutput::PPMOutput(SP_mesh mesh)
  : d_mesh(mesh)
{
  Require(d_mesh);
}

//---------------------------------------------------------------------------//
bool PPMOutput::initialize(const std::string filename, size_t n)
{
  Require(filename != "");
  Require(n > 0);

  d_filename = filename;

  // physical extent
  double X = d_mesh->total_width_x();
  double Y = d_mesh->total_width_y();

  // minimum pixel size
  d_dxyx = detran_utilities::vec_min(d_mesh->dx());
  d_dxyx = std::min(d_dxyx, detran_utilities::vec_min(d_mesh->dy()));

  // resolution
  d_nx = n*std::ceil(X/d_dxyx);
  d_ny = n*std::ceil(Y/d_dxyx);
  d_x = detran_utilities::linspace(0, d_mesh->total_width_x(), d_nx);
  d_y = detran_utilities::linspace(0, d_mesh->total_width_y(), d_ny);
  d_cells.resize(d_nx * d_ny);
  for (size_t j = 0; j < d_ny; ++j)
  {
    for (size_t i = 0; i < d_nx; ++i)
    {
      d_cells[i + j * d_ny] =
        d_mesh->find_cell(detran_utilities::Point(d_x[i], d_y[j]));
      if (d_cells[i + j * d_ny] < 0.0)
      {
        std::cout << " missing cell for x = " << d_x[i] << " y = " << d_y[j] << std::endl;
      }
    }
  }

  // create plotter
  d_plotter = new PPMPlotter;

  return true;
}

//---------------------------------------------------------------------------//
bool PPMOutput::write_mesh_map(const std::string &key)
{
  if (!d_mesh->mesh_map_exists(key))
  {
    std::cout << "Mesh map " + key + " doesn't exist; not writing mesh map"
              << std::endl;
    return false;
  }
  d_plotter->initialize(d_nx, d_ny, d_filename+"_"+key+".ppm");
  vec_dbl data(d_cells.size());
  for (size_t i = 0; i < d_cells.size(); ++i)
  {
    Assert(d_cells[i] >= 0.0);
    data[i] = d_mesh->mesh_map(key)[d_cells[i]];
  }
  d_plotter->set_pixels<double>(data);
  return d_plotter->write();
}

//---------------------------------------------------------------------------//
bool PPMOutput::write_scalar_flux(SP_state state)
{
  for (size_t g = 0; g < state->number_groups(); ++g)
  {
    std::cout << " g=" << g << std::endl;
    std::stringstream g_;
    g_ << g;
    d_plotter->initialize(d_nx, d_ny, d_filename+"_scalar_flux_"+g_.str()+".ppm");
    vec_dbl data(d_cells.size());
    for (size_t i = 0; i < d_cells.size(); ++i)
    {
      //Assert(d_cells[i] >= 0.0);
      data[i] = state->phi(g)[d_cells[i]];
    }
    d_plotter->set_pixels<double>(data);
    d_plotter->write();
  }
  return true;
}

} // end namespace detran_ioutils



