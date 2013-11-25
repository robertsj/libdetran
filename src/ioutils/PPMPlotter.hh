//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PPMPlotter.hh
 *  @brief PPMPlotter class definition
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_ioutils_PPMPLOTTER_HH_
#define detran_ioutils_PPMPLOTTER_HH_

#include "ioutils/ioutils_export.hh"
#include "ioutils/ColorMap.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <string>

namespace detran_ioutils
{

/**
 *  @class PPMPlotter
 *  @brief Produces 2-D plots in the PPM format
 *  @todo  Add capability to do xz and yz planes for 3-D
 */
class IOUTILS_EXPORT PPMPlotter
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<PPMPlotter>      SP_ppmplotter;
  typedef detran_utilities::size_t              size_t;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef detran_utilities::vec_int             vec_int;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //--------------------------------------------------------------------------//

  PPMPlotter();

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Initialize an image of nx by ny pixels
  void initialize(const size_t       nx,
                  const size_t       ny,
                  const std::string &name = "plot.ppm",
                  const size_t       scheme = ColorMap::DEFAULT);

  /**
   *  @brief Set the value for pixel (i, j)
   *
   *  The user sets the real value.  The data is scaled appropriately
   *  between the minimum and maximum to maximize use of the available
   *  color range (8 bit).  Note, the value must be *positive*.
   */
  void set_pixel(const size_t i, const size_t j, const double v);

  /// Write all pixels at once.
  template <class T>
  void set_pixels(const std::vector<T>& v)
  {
    Require(v.size() == d_image.size());
    for (size_t i = 0; i < d_nx * d_ny; ++i) d_image[i] = v[i];
  }

  /// Write the file
  bool write();

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Horizontal resolution
  size_t d_nx;
  /// Vertical resolution
  size_t d_ny;
  /// Image data as 1-D array
  vec_dbl d_image;
  /// Name
  std::string d_name;
  /// Color scheme
  size_t d_scheme;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

IOUTILS_TEMPLATE_EXPORT(detran_utilities::SP<PPMPlotter>)

} // end namespace detran_ioutils

#endif // detran_ioutils_PPMPLOTTER_HH_

//----------------------------------------------------------------------------//
//              end of file PPMPlotter.hh
//----------------------------------------------------------------------------//
