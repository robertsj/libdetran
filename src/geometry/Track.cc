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
Track::Track(const Point  &r0, const Point  &r1, double w)
  : d_enter(r0)
  , d_exit(r1)
  , d_spatial_weight(w)
  , d_identifier(d_number_tracks++)
  , d_reversed(false)
  , d_reversed_identifier(d_identifier)
  , d_segments(d_segments_)
{

}

//----------------------------------------------------------------------------//
Track::SP_track Track::reverse()
{
  return SP_track(new Track(*this, true));
}


//----------------------------------------------------------------------------//
Track::Track(Track &t0, bool reversed)
  : d_enter(t0.enter())
  , d_exit(t0.exit())
  , d_spatial_weight(t0.spatial_weight())
  , d_identifier(d_number_tracks++)
  , d_reversed(reversed)
  , d_reversed_identifier(t0.identifier())
  , d_segments(t0.segments())
{

}

//----------------------------------------------------------------------------//
const Point& Track::enter() const
{
  if (d_reversed)
    return d_exit;
  return d_enter;
}

//----------------------------------------------------------------------------//
const Point&  Track::exit() const
{
  if (d_reversed)
    return d_enter;
  return d_exit;
}

//----------------------------------------------------------------------------//
void Track::add_segment(const Segment &s)
{
  Require(d_reversed == false);
  d_segments.push_back(s);
}

//----------------------------------------------------------------------------//
int Track::number_segments() const
{
  return d_segments.size();
}

//----------------------------------------------------------------------------//
Track::iterator Track::begin()
{
  return iterator(&d_segments, not d_reversed).begin();
}

//----------------------------------------------------------------------------//
Track::iterator Track::end()
{
  return iterator(&d_segments, not d_reversed).end();
}

//----------------------------------------------------------------------------//
const int& Track::identifier() const
{
  return d_identifier;
}

//----------------------------------------------------------------------------//
double Track::spatial_weight() const
{
  return d_spatial_weight;
}

//----------------------------------------------------------------------------//
bool Track::operator< (const Track &t2) const
{
  if (enter() <= t2.enter())
   {
     return true;
   }
   else if (enter() == t2.enter())
   {
     return exit() < t2.exit();
   }
  return false;
}

//----------------------------------------------------------------------------//
bool Track::operator== (const Track &t2) const
{
  if ((enter() == t2.enter()) and ((exit() == t2.exit())))
    return true;
   return false;
}

//----------------------------------------------------------------------------//
std::ostream& operator<< (std::ostream &out, Track &t)
{
  std::ios::fmtflags f(out.flags());
  out.setf(std::ios::fixed, std::ios::floatfield);
  out.setf(std::ios::showpoint);
  out << "    id = " << t.identifier() << std::endl;
  out << " enter = " << t.enter() << std::endl
      << "  exit = "  << t.exit() << std::endl;
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
