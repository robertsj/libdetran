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
Region::Region()
{
  /* ... */
}

//----------------------------------------------------------------------------//
Region::SP_region Region::Create()
{
  SP_region p(new Region());
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

} // end namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of Region.cc
//----------------------------------------------------------------------------//
