//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CSG.cc
 *  @brief Member definitions of CSG-related classes
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "CSG.hh"

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
CSG_Primitive::intersections(c_Point &r, c_Point &d, c_double t_max)
{
//  std::cout << " prim r = " << r << std::endl;
  CSG_Primitive::vec_point points =  d_surface->intersections(r, d);
//  for (int i = 0; i < points.size(); ++i)
//    std::cout << " prim intersection " << i << " " << points[i] << std::endl;
  return points;
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
CSG_Operator::intersections(c_Point &r, c_Point &d, c_double t_max)
{
//  std::cout << "op_int r = " << r  << " d = " << d << std::endl;

  // get intersections with left node and eliminate those not contained
  vec_point pointsL = d_L->intersections(r, d, t_max);

//  for (int i = 0; i < pointsL.size(); ++i)
//    std::cout << " opL intersection " << i << " " << pointsL[i] << std::endl;

//  for (vec_point::iterator p = pointsL.begin(); p != pointsL.end(); ++p)
//    if (!contains(*p)) p = pointsL.erase(p);

//  for (int i = 0; i < pointsL.size(); ++i)
//    std::cout << " opL intersection2 " << i << " " << pointsL[i] << std::endl;

  // do the same for the right
  vec_point pointsR = d_R->intersections(r, d, t_max);

//  for (int i = 0; i < pointsR.size(); ++i)
//    std::cout << " opR intersection " << i << " " << pointsR[i] << std::endl;


//  for (vec_point::iterator p = pointsR.begin(); p != pointsR.end(); ++p)
//    if (!contains(*p)) p = pointsR.erase(p);

//  for (int i = 0; i < pointsR.size(); ++i)
//    std::cout << " opR intersection 2 " << i << " " << pointsR[i] << std::endl;


  // concatenate the results and sort by increasing distance from origin r
  pointsL.insert(pointsL.end(), pointsR.begin(), pointsR.end());
  std::sort(pointsL.begin(), pointsL.end(), IntersectionPointCompare(r));

//  for (int i = 0; i < pointsL.size(); ++i)
//    std::cout << " op final intersection " << i << " " << pointsL[i] << std::endl;

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
  return (d_L->contains(r) && !d_R->contains(r)) ||
         (d_R->contains(r) && !d_L->contains(r));
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
CSG_Translation::intersections(c_Point &r, c_Point &d, c_double t_max)
{
  return d_node->intersections(r - d_translation, d, t_max);
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
