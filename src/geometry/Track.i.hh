//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  Track.i.hh
 *  @brief Track inline member definitions
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_geometry_TRACK_I_HH_
#define detran_geometry_TRACK_I_HH_

namespace detran_geometry
{

//----------------------------------------------------------------------------//
inline Track::Point Track::enter() const
{
  return d_enter;
}

//----------------------------------------------------------------------------//
inline Track::Point Track::exit() const
{
  return d_exit;
}

//----------------------------------------------------------------------------//
inline Track::size_t Track::number_segments() const
{
  return d_segments.size();
}

//----------------------------------------------------------------------------//
inline Segment& Track::segment(const size_t i)
{
  Require(i < d_segments.size());
  return d_segments[i];
}

//----------------------------------------------------------------------------//
inline const Segment& Track::segment(const size_t i) const
{
  Require(i < d_segments.size());
  return d_segments[i];
}

//----------------------------------------------------------------------------//
inline Track::iterator Track::begin()
{
  return d_segments.begin();
}

//----------------------------------------------------------------------------//
inline Track::riterator Track::rbegin()
{
  return d_segments.rbegin();
}

//----------------------------------------------------------------------------//
inline double Track::cos_phi() const
{
  return d_cos_phi;
}

//----------------------------------------------------------------------------//
inline double Track::sin_phi() const
{
  return d_sin_phi;
}

//----------------------------------------------------------------------------//
inline double Track::width() const
{
  return d_width;
}

} // end namespace detran_geometry

#endif /* detran_geometry_TRACK_I_HH_ */

//----------------------------------------------------------------------------//
//              end of file Track.i.hh
//----------------------------------------------------------------------------//
