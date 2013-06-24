//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PPMPlotter.cc
 *  @brief PPMPlotter member definitions
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//----------------------------------------------------------------------------//

#include "ioutils/PPMPlotter.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/Random.hh"
#include <fstream>
#include <cstdio>

namespace detran_ioutils
{

//----------------------------------------------------------------------------//
PPMPlotter::PPMPlotter()
  : d_nx(0)
  , d_ny(0)
  , d_scheme(ColorMap::DEFAULT)
{
  /* ... */
}

//----------------------------------------------------------------------------//
void PPMPlotter::
initialize(const size_t        nx,
           const size_t        ny,
           const std::string  &name,
           const size_t        scheme)
{
  Require(nx > 0);
  Require(ny > 0);
  Require(name.length() > 0);
  Require(scheme < ColorMap::END_COLOR_MAPS);
  d_nx = nx;
  d_ny = ny;
  d_name = name;
  d_image.resize(nx * ny, -1.0);
  d_scheme = scheme;
}

//----------------------------------------------------------------------------//
void PPMPlotter::set_pixel(const size_t i, const size_t j, const double v)
{
  d_image[i + j * d_nx] = v;
}

//----------------------------------------------------------------------------//
bool PPMPlotter::write()
{
  // open the file
  std::ofstream file;
  file.open(d_name.c_str(), std::ios::binary);

  // write the header
  file << "P6\n" << d_nx << " " << d_ny << "\n" << 255 << "\n";

  // get the colors
  ColorMap::vec_rgb colors = ColorMap::color(d_scheme, d_image);

  // write each pixel
  for (size_t i = 0; i < d_nx * d_ny; ++i)
  {
    file << colors[i].r << colors[i].g << colors[i].b;
  }
  file.close();
  return true;
}

} // end namespace detran_ioutils

//----------------------------------------------------------------------------//
//              end of file PPMPlotter.cc
//----------------------------------------------------------------------------//
