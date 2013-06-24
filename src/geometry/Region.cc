//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Region.cc
 *  @brief Region member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Region.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Region::Region(const size_t mat, const Point &b_min, const Point &b_max)
  : d_bounds(b_min, b_max)
  , d_have_bound(false)
  , d_have_bound_z(false)
{
  if (distance(d_bounds[0], d_bounds[1]) > 0.0)
    d_have_bound = true;
  if (d_have_bound)
  {
    Require(d_bounds[0] <= d_bounds[1]);
    if (d_bounds[1].z() > d_bounds[0].z())
      d_have_bound_z = true;
  }
  add_attribute("MATERIAL", mat);
}

//----------------------------------------------------------------------------//
Region::SP_region
Region::Create(const size_t mat, const Point &b_min, const Point &b_max)
{
  SP_region p(new Region(mat, b_min, b_max));
  return p;
}

//----------------------------------------------------------------------------//
void Region::append(SP_region r, const size_t op)
{
  Require(r);
  Require(r->top_node());
  append_node(d_node, r->top_node(), op);
}

//----------------------------------------------------------------------------//
void Region::append(SP_node n, const size_t op)
{
  Require(n);
  append_node(d_node, n, op);
}

//----------------------------------------------------------------------------//
void Region::append(SP_surface s, bool sense)
{
  Require(s);
  SP_node n(new CSG_Primitive(s, sense));
  // Surfaces are only added via intersection, as sense contains enough
  append_node(d_node, n, INTERSECTION);
}

//----------------------------------------------------------------------------//
bool Region::contains(const Point &r)
{
  Require(d_node);

  // If we have a bounding box, eliminate any points not in the box
  if (d_have_bound)
    if (r > d_bounds[0] && r < d_bounds[1])
      return d_node->contains(r);
  return d_node->contains(r);
}

//----------------------------------------------------------------------------//
bool Region::intersects_bounding_box(const Ray &r, const double max_length)
{
  // This is nearly verbatim to the code given in:
  //
  // Williams et al. "An Efficient and Robust Ray-Box Intersection
  //   Algorithm", Journal of Graphics, GPU, and Game Tools, 10 (2005)
  //
  double tmin, tmax, tymin, tymax, tzmin, tzmax;
  tmin  = (d_bounds[    r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
  tmax  = (d_bounds[1 - r.sign[0]].x() - r.origin.x()) * r.inv_direction.x();
  tymin = (d_bounds[    r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
  tymax = (d_bounds[1 - r.sign[1]].y() - r.origin.y()) * r.inv_direction.y();
  if ((tmin > tymax) || (tymin > tmax))
    return false;
  if (tymin > tmin)
    tmin = tymin;
  if (tymax < tmax)
    tmax = tymax;
  if (d_have_bound_z)
  {
    tzmin = (d_bounds[    r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
    tzmax = (d_bounds[1 - r.sign[2]].z() - r.origin.z()) * r.inv_direction.z();
    if ((tmin > tzmax) || (tzmin > tmax))
      return false;
    if (tzmin > tmin)
      tmin = tzmin;
    if (tzmax < tmax)
      tmax = tzmax;
  }
  return ((tmin < max_length) && (tmax > 0.0));
}

//----------------------------------------------------------------------------//
void Region::append_node(SP_node current, SP_node addition, const size_t op)
{
  Require(addition);
  if (!current)
  {
    // This is the first node added.
    d_node = addition;
    return;
  }
  Require(op < END_NODE_OPERATORS);
  if (op == UNION)
    d_node = new CSG_Union(current, addition);
  else if (op == INTERSECTION)
    d_node = new CSG_Intersection(current, addition);
  else
    d_node = new CSG_Difference(current, addition);
}

//----------------------------------------------------------------------------//
void Region::add_attribute(const std::string &key, const int value)
{
  Require(!key.empty());
  // Erase the value associated with the key if it exists.
  attributes_t::iterator iter;
  iter = d_attributes.find(key);
  if (iter != d_attributes.end())
    d_attributes.erase(iter);
  // Add the new value.
  d_attributes[key] = value;
}

//----------------------------------------------------------------------------//
bool Region::attribute_exists(const std::string &key) const
{
  Require(!key.empty());
  attributes_t::const_iterator iter;
  iter = d_attributes.find(key);
  if (iter != d_attributes.end())
    return true;
  else
    return false;
}

//----------------------------------------------------------------------------//
int Region::attribute(const std::string &key) const
{
  Require(!key.empty());
  attributes_t::const_iterator iter;
  iter = d_attributes.find(key);
  Ensurev(iter != d_attributes.end(), "Region attribute not found.");
  return iter->second;
}

//----------------------------------------------------------------------------//
Point Region::bound_min() const
{
  return d_bounds[0];
}

//----------------------------------------------------------------------------//
Point Region::bound_max() const
{
  return d_bounds[1];
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Region.cc
//----------------------------------------------------------------------------//
