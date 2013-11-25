//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ColorMap.hh
 *  @brief ColorMap class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_ioutils_COLORMAP_HH_
#define detran_ioutils_COLORMAP_HH_

#include "geometry/Point.hh"
#include "utilities/Definitions.hh"
#include "utilities/MathUtilities.hh"
#include "utilities/Random.hh"

namespace detran_ioutils
{

/**
 *  @class ColorMap
 *  @brief Provides basic value-to-color mapping for simple plotting
 */
class ColorMap
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//

  /// Available color maps
  enum COLOR_MAPS
  {
    DEFAULT,        /// blue-greenish-red
    COOLWARM,       /// "diverging" blue-white-red
    JET,            /// close to the matlab jet
    HOT,            /// close to the matlab hot
    COOL,           /// close to the matlab cool
    RANDOM,         /// random colors for each unique value
    END_COLOR_MAPS
  };

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_dbl         vec_dbl;
  struct rgb_t
  {
    unsigned char r, g, b;
  };
  typedef std::vector<ColorMap::rgb_t>      vec_rgb;
  typedef detran_geometry::Point            Point;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Return a vector of colors representing the values
   *  @param  map     Color map identifier
   *  @param  values  Values to represent on the color map
   */
  static vec_rgb color(const size_t map, const vec_dbl &values);

  static ColorMap::rgb_t hex_to_rgb(const int hex);

private:

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

  /// Color maps defined as functions of value
  static vec_rgb colors_function(const size_t map, const vec_dbl &values);
  /// Random color for unique values
  static vec_rgb colors_random(const vec_dbl &values);
  /// Color maps defined as interpolation between values
  static vec_rgb colors_interp(const size_t map, const vec_dbl &values);

public:

  /// Basic 16 HTML colors + a few extras
  enum COLORS
  {
    azure =       0x007FFF,
    black =       0x000000,
    blue =        0x0000FF,
    cyan =        0x00FFFF,
    darkblue =    0x00008B,
    darkred =     0x8B0000,
    fuchsia =     0xFF00FF,
    gray =        0x808080,
    green =       0x008000,
    lightgreen =  0x90EE90,
    lime =        0x00FF00,
    maroon =      0x800000,
    navy =        0x000080,
    olive =       0x808000,
    orange =      0xFFA500,
    purple =      0x800080,
    red =         0xFF0000,
    silver =      0xC0C0C0,
    teal =        0x008080,
    white =       0xFFFFFF,
    yellow =      0xFFFF00
  };

};

//----------------------------------------------------------------------------//
struct gray_to_color
{
  bool operator==(const double v) const
  {
    return v == value;
  }
  double            value;
  ColorMap::rgb_t   color;
};

} // end namespace detran_ioutils

#endif /* detran_ioutils_COLORMAP_HH_ */

//----------------------------------------------------------------------------//
//              end of file ColorMap.cc
//----------------------------------------------------------------------------//
