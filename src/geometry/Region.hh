//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Region.hh
 *  @brief Region class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_REGION_HH_
#define detran_geometry_REGION_HH_

#include "CSG.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"
#include <map>

namespace detran_geometry
{

/**
 *  @class Region
 *  @brief Represents a constant material region described via CSG
 *
 *  A solid is comprised of a combination of primitives, in this case
 *  surfaces, via the union, intersection, and difference operations.
 *  Regions can be created directly from surfaces, nodes, and other regions.
 */
class Region
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef detran_utilities::SP<Region>      SP_region;
  typedef std::vector<SP_region>            vec_region;
  typedef Surface::SP_surface               SP_surface;
  typedef Surface::vec_surface              vec_surface;
  typedef CSG_Node::SP_node                 SP_node;
  typedef detran_utilities::size_t          size_t;
  typedef detran_utilities::vec_dbl         vec_dbl;
  typedef detran_utilities::vec_size_t      vec_size_t;
  typedef std::map<std::string, int>        attributes_t;

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// Constructor with optional bounding box to expedite tracking
  Region(const size_t mat, const vec_dbl &bbox = vec_dbl(0));
  /// SP constructor
  static SP_region Create(const size_t mat, const vec_dbl &bbox = vec_dbl(0));
  /// Add a region
  void append(SP_region region, const size_t op);
  /// Add a new node
  void append(SP_node node, const size_t op);
  /// Add a new surface node
  void append(SP_surface surface, bool sense);
  /// Does the region contain the point?
  bool contains(const Point &r);
  /// Return the top node
  SP_node top_node() {return d_node;}
  /// Add an attribute to this region
  void add_attribute(const std::string &key, const int value);
  /// Get a region attribute
  int attribute(const std::string &key) const;
  /// Check for an attribute
  bool attribute_exists(const std::string &key) const;
  /// Get the bounding box
  vec_dbl bounding_box() const;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// CSG node
  SP_node d_node;

  /**
   *  Bounding box that can be used to expedite tracking and plotting.  This
   *  must be *larger* than any boundary of the region, i.e. *no* overlap.
   *  Unfortunately, the onus is on the user to give a correct bounding box if
   *  given at all.
   */
  vec_dbl d_bounding_box;

  /// String-keyed integer attributes for the region
  attributes_t d_attributes;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  void append_node(SP_node current, SP_node addition, const size_t op);

};

} // end namespace detran_geometry

#endif /* detran_geometry_REGION_HH_ */

//----------------------------------------------------------------------------//
//              end of Region.hh
//----------------------------------------------------------------------------//
