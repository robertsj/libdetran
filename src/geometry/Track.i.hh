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
inline Point Track::enter() const
{
  return d_enter;
}

//----------------------------------------------------------------------------//
inline Point Track::exit() const
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
inline Track::iterator Track::begin(bool forward)
{
  return iterator(&d_segments, forward);
}

//----------------------------------------------------------------------------//
inline Track::iterator Track::end(bool forward)
{
  return iterator(&d_segments, forward) + d_segments.size();
}

//----------------------------------------------------------------------------//
inline double Track::mu() const
{
  return d_mu;
}

//----------------------------------------------------------------------------//
inline double Track::eta() const
{
  return d_eta;
}

//----------------------------------------------------------------------------//
inline double Track::xi() const
{
  return d_xi;
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
