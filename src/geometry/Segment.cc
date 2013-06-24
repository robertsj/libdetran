//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Segment.cc
 *  @brief Segment class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Segment.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Segment::Segment(const size_t r, const double l)
  : d_region(r)
  , d_length(l)
{
  Require(l >= 0.0);
}

//----------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream &out, const Segment &s)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "(region = " << s.region() << ", length = " << s.length() << ")";
  out.flags(f);
  return out;
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Segment.cc
//----------------------------------------------------------------------------//



