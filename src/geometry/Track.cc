//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Track.cc
 *  @brief Track member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "Track.hh"

namespace detran_geometry
{

//----------------------------------------------------------------------------//
Track::size_t Track::d_number_tracks = 0;

//----------------------------------------------------------------------------//
Track::Track(const Point &r0,
             const Point &r1,
             const double w)
  : d_enter(r0)
  , d_exit(r1)
  , d_width(w)
  , d_identifier(d_number_tracks++)
{
  double d = distance(r0, r1);
  Point p = r1 - r0;
  d_mu  = p.x() / d;
  d_eta = p.y() / d;
  d_xi  = p.z() / d;
}

//----------------------------------------------------------------------------//
Track::SP_track Track::Create(const Point  &r0,
                              const Point  &r1,
                              const double  w)
{
  SP_track p(new Track(r0, r1, w));
  return p;
}

//----------------------------------------------------------------------------//
void Track::add_segment(const Segment &s)
{
  d_segments.push_back(s);
}

//----------------------------------------------------------------------------//
Track::size_t Track::identifier() const
{
  return d_identifier;
}

//----------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream &out, Track &t)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << " enter = " << t.enter()
      << " exit = "  << t.exit() << std::endl;
  for (int s = 0; s < t.number_segments(); s++)
  {
    out << "          segment = " << s << " " << t.segment(s) << std::endl;
  }
  out.flags(f);
  return out;
}

} // namespace detran_geometry

//----------------------------------------------------------------------------//
//              end of file Track.cc
//----------------------------------------------------------------------------//
