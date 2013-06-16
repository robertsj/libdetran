//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  Geometry.hh
 *  @brief Geometry class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_GEOMETRY_HH_
#define detran_GEOMETRY_HH_

#include "geometry/Region.hh"

namespace detran_geometry
{

/**
 *  @class Geometry
 *  @brief Represents a complete geometry comprised of @ref Region objects
 */
class GEOMETRY_EXPORT Geometry
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Geometry>    SP_geometry;
  typedef Region::SP_region                 SP_region;
  typedef std::vector<SP_region>            vec_region;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_size_t      vec_size_t;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Constructor
  Geometry(const double x, const double y, const double z);
  /// SP constructor
  static SP_geometry Create(const double x, const double y, const double z);
  /// Add a region
  void add_region(SP_region r);
  /// Get region
  SP_region region(const size_t r);
  /// Find region containing cell. Returns index >= 0 or -1 if none found
  int find(const Point &r);
  /// Return number of regions in the system
  size_t number_regions() const;
  /// Return material index for a region
  size_t material_index(const size_t r) const;
  //@{
  /// Get geometry bounding box coordinates
  double width_x() const {return d_x;}
  double width_y() const {return d_y;}
  double width_z() const {return d_z;}
  //@}
  size_t dimension() const {return 2;}

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  //@{
  ///  Bounding box coordinates, where the origin is always at [0, 0, 0]
  double d_x;
  double d_y;
  double d_z;
  //@}
  /// Regions
  vec_region d_regions;
  /// Map specifying materials in each region
  vec_size_t d_material_map;

};

} // end namespace detran_geometry

#endif /* detran_GEOMETRY_HH_ */

//----------------------------------------------------------------------------//
//              end of Geometry.hh
//----------------------------------------------------------------------------//
