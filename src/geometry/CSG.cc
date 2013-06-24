//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CSG.cc
 *  @brief Member definitions of CSG-related classes
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CSG.hh"
#include <algorithm>

namespace detran_geometry
{

//----------------------------------------------------------------------------//
// PRIMITIVE
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Primitive::CSG_Primitive(SP_surface surface, bool sense)
  : d_surface(surface)
  , d_sense(sense)
{
  Require(surface);
}

//----------------------------------------------------------------------------//
bool CSG_Primitive::contains(c_Point &r) const
{
  return d_surface->sense(r) == d_sense;
}

//----------------------------------------------------------------------------//
CSG_Primitive::vec_point
CSG_Primitive::intersections(c_Ray &r, c_dbl t_max)
{
  return d_surface->intersections(r, t_max);
}

//----------------------------------------------------------------------------//
// OPERATOR
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Operator::CSG_Operator(SP_node L, SP_node R)
  : d_L(L), d_R(R)
{
  Require(d_L);
  Require(d_R);
  Require(d_R != d_L);
}

//----------------------------------------------------------------------------//
CSG_Operator::vec_point
CSG_Operator::intersections(c_Ray &r, c_dbl t_max)
{
  // get intersections with left and right nodes
  vec_point pointsL = d_L->intersections(r, t_max);
  vec_point pointsR = d_R->intersections(r, t_max);

  // concatenate the results and sort by increasing distance from origin r
  pointsL.insert(pointsL.end(), pointsR.begin(), pointsR.end());
  std::sort(pointsL.begin(), pointsL.end(), IntersectionPointCompare(r.origin));

  return pointsL;
}

//----------------------------------------------------------------------------//
// UNION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Union::CSG_Union(SP_node L, SP_node R)
  : CSG_Operator(L, R)
{
  /* ... */
}

//----------------------------------------------------------------------------//
bool CSG_Union::contains(const Point &r) const
{
  return d_L->contains(r) || d_R->contains(r);
}

//----------------------------------------------------------------------------//
// INTERSECTION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Intersection::CSG_Intersection(SP_node L, SP_node R)
  : CSG_Operator(L, R)
{
  /* ... */
}

//----------------------------------------------------------------------------//
bool CSG_Intersection::contains(const Point &r) const
{
  return d_L->contains(r) && d_R->contains(r);
}

//----------------------------------------------------------------------------//
// DIFFERENCE
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Difference::CSG_Difference(SP_node L, SP_node R)
  : CSG_Operator(L, R)
{
  /* ... */
}

//----------------------------------------------------------------------------//
bool CSG_Difference::contains(const Point &r) const
{
  // Difference consists of those points in L but not in R.  It's equivalent
  // to (L .intersection. (.not. R)), but we don't have a .not. operator.
  return d_L->contains(r) && !d_R->contains(r);
}

//----------------------------------------------------------------------------//
// TRANSLATION
//----------------------------------------------------------------------------//

//----------------------------------------------------------------------------//
CSG_Translation::CSG_Translation(SP_node node, c_Point &translation)
  : d_node(node)
  , d_translation(translation)
{
  Require(node);
}

//----------------------------------------------------------------------------//
CSG_Translation::vec_point
CSG_Translation::intersections(c_Ray &r, c_dbl t_max)
{
  Ray r_tran = Ray(r.origin-d_translation, r.direction);
  vec_point points = d_node->intersections(r_tran, t_max);
  for (size_t i = 0; i < points.size(); ++i)
    points[i] = points[i] + d_translation;
  return points;
}

//----------------------------------------------------------------------------//
bool CSG_Translation::contains(const Point &r) const
{
  return d_node->contains(r - d_translation);
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file CSG.cc
//----------------------------------------------------------------------------//
