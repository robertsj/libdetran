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
Region::Region(const size_t mat, const vec_dbl &bbox)
  : d_bounding_box(bbox)
{
  if (d_bounding_box.size())
  {
    Require(d_bounding_box.size() == 6);
  }
  add_attribute("MATERIAL", mat);
}

//----------------------------------------------------------------------------//
Region::SP_region Region::Create(const size_t mat, const vec_dbl &bbox)
{
  SP_region p(new Region(mat, bbox));
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
  if (d_bounding_box.size())
  {
    if (r.x() > d_bounding_box[0] && r.x() < d_bounding_box[1] &&
        r.y() > d_bounding_box[2] && r.y() < d_bounding_box[3] &&
        r.z() > d_bounding_box[4] && r.z() < d_bounding_box[5])
    {
      return d_node->contains(r);
    }
  }
  return d_node->contains(r);
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
Region::vec_dbl Region::bounding_box() const
{
  return d_bounding_box;
}

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Region.cc
//----------------------------------------------------------------------------//
