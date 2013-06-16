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

  // Get the colors
  ColorMap::vec_rgb colors = ColorMap::color(d_scheme, d_image);

  // scale the data
//  double min = detran_utilities::vec_min(d_image);
//  double max = detran_utilities::vec_max(d_image);
//  detran_utilities::vec_plus_a(d_image, -min);
//  if (max-min > 0) detran_utilities::vec_scale(d_image, 1.0/(max-min));

  // write each pixel
  for (size_t i = 0; i < d_nx * d_ny; ++i)
  {
    // scale data to [0, 1], and then define r, g, and b.
    // this scale should give fair cold-to-hot.  Default
    // is white, which indicates the pixel is "outside"
    // whatever is being plotted.
//    unsigned char r = 255;
//    unsigned char g = 255;
//    unsigned char b = 255;
//    if (d_image[i] >= 0.0)
//    {
//      double v = (d_image[i]-min)/(max-min);
//      r = 255. * (2.0*v - v*v);
//      g = 255. * (4.0*v - 4.0*v*v);
//      b = 255. * (1.0-v*v);
//    }
//    //printf("%4i %4i %4i %4i %8.3f \n", i, r, g, b, d_image[i]);
//    file << r << g << b;
    file << colors[i].r << colors[i].g << colors[i].b;
  }
  file.close();
  return true;
}

//PPMPlotter::rgb_t PPMPlotter::color(const size_t scheme, const double v)
//{
//  rgb_t rgb = {0, 0, 0};
//  if (scheme == COLDHOT)
//  {
//    rgb[0] = 255. * (2.0*v - v*v);
//    rgb[1] = 255. * (4.0*v - 4.0*v*v);
//    rgb[2] = 255. * (1.0-v*v);
//  }
//  else
//  {
//
//  }
//  return rgb;
//}

} // end namespace detran_ioutils

//----------------------------------------------------------------------------//
//              end of file PPMPlotter.cc
//----------------------------------------------------------------------------//
