//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PPMPlotter.hh
 *  @brief  PPMPlotter 
 *  @author Jeremy Roberts
 *  @date   Mar 12, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_ioutils_PPMPLOTTER_HH_
#define detran_ioutils_PPMPLOTTER_HH_

#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"

#include <string>

namespace detran_ioutils
{

/**
 *  @class PPMPlotter
 *  @brief Produces 2-D plots in the PPM format
 *
 *  This is serves as a fast way to plot slices of the
 *  geometry.
 */
class PPMPlotter
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::size_t        size_t;
  typedef detran_utilities::vec_dbl       vec_dbl;
  typedef detran_utilities::vec_int       vec_int;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //-------------------------------------------------------------------------//

  PPMPlotter();

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Initialize an image of nx by ny pixels
  void initialize(const size_t nx,
                  const size_t ny,
                  std::string name = "plot.ppm");

  /**
   *  @brief Set the value for pixel (i, j)
   *
   *  The user sets the real value.  The data is scaled appropriately
   *  between the minimum and maximum to maximize use of the available
   *  color range (8 bit).  Note, the value must be *positive*.
   */
  void set_pixel(const size_t i, const size_t j, const double v);

  /// Write the file
  bool write();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Horizontal resolution
  size_t d_nx;
  /// Vertical resolution
  size_t d_ny;
  /// Image data as 1-D array
  vec_dbl d_image;
  /// Name
  std::string d_name;

};

} // end namespace detran_ioutils

#endif // detran_ioutils_PPMPLOTTER_HH_

//---------------------------------------------------------------------------//
//              end of file PPMPlotter.hh
//---------------------------------------------------------------------------//
