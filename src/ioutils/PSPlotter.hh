//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  PSPlotter.hh
 *  @brief PSPlotter class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//


#ifndef detran_ioutils_PSPLOTTER_HH_
#define detran_ioutils_PSPLOTTER_HH_

#include "ioutils/ioutils_export.hh"
#include "ioutils/ColorMap.hh"
#include "geometry/Point.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <string>
#include <fstream>

namespace detran_ioutils
{

/**
 *  @class PSPlotter
 *  @brief Plot lines in PS format
 */
class IOUTILS_EXPORT PSPlotter
{

public:

  //--------------------------------------------------------------------------//
  // ENUMERATIONS
  //--------------------------------------------------------------------------//


  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::size_t              size_t;
  typedef detran_utilities::vec_dbl             vec_dbl;
  typedef detran_utilities::vec_int             vec_int;
  typedef detran_geometry::Point                Point;


  //--------------------------------------------------------------------------//
  // CONSTRUCTOR AND DESTRUCTOR
  //--------------------------------------------------------------------------//

//  %!PS-Adobe-3.0 EPSF-3.0
//  %%Document-Fonts: Times-Roman
//  %%Title: circle.eps
//  %%Creator: PS_Write.F
//  %%CreationDate: 02-Aug-99
//  %%Pages: 1
//  %%BoundingBox:   36   36  576  756


  PSPlotter(std::string name, vec_dbl bbox)
    : d_name(name)
    , d_bounding_box(bbox)
  {
    d_file.open(d_name.c_str());
    d_file << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl
           << "%%Creator: Detran, a deterministic transport code" << std::endl
           << "%%BoundingBox: "
           << d_bounding_box[0] << " " << d_bounding_box[1] << " "
           << d_bounding_box[2] << " " << d_bounding_box[3] << std::endl;
  }

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  void draw_point(const Point &P0)
  {

  }
  void draw_line(const Point &P0, const Point &P1)
  {
    d_file << "newpath" << std::endl
           << 100.0*P0.x() << " " << 100.0*P0.y() << " moveto " << std::endl
           << 100.0*P1.x() << " " << 100.0*P1.y() << " lineto " << std::endl
           << "stroke" << std::endl;
  }
  void finalize()
  {
    d_file.close();
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Name
  std::string d_name;
  /// Bounding box
  vec_dbl d_bounding_box;
  /// Filestream
  std::ofstream d_file;

  //--------------------------------------------------------------------------//
  // IMPLEMENTATION
  //--------------------------------------------------------------------------//

};

} // end namespace detran_ioutils

#endif /* detran_ioutils_PSPLOTTER_HH_ */
