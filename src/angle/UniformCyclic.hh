//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  UniformCyclic.hh
 *  @brief UniformCyclic class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_UNIFORM_HH_
#define detran_angle_UNIFORM_HH_

#include "ProductQuadrature.hh"

namespace detran_angle
{

/**
 *  @class UniformCyclic
 *  @brief Azimuthal quadrature set for cyclic tracking in 2-D MOC
 *
 *  Cyclic tracking is not explicitly implemented in Detran, but by
 *  adjusting azimuths so that tracking is cyclic, any interpolation at
 *  the boundary becomes essentially exact for the case of even track
 *  spacing.  This is limited to the case of uniformly-spaced tracks in
 *  square cells.  The client must know a priori how many tracks are
 *  created for each azimuth in a quadrant, and hence, this is mostly
 *  a class for development.
 *
 */
class ANGLE_EXPORT UniformCyclic
{

public:

  //--------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //--------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param dim                  Problem dimension (only 2 supported for now)
   *  \param num_azimuths_octant  Number of azimuths per octant
   *  \param num_space            Number of tracks per azimuth
   *  \param num_polar            Number of polar angles in half space
   *  \param polar                Polar quadrature string identifier
   */
  UniformCyclic(size_t dim,
          size_t num_azimuths_octant,
          size_t num_space,
          size_t num_polar,
          std::string polar);

  /// SP Constructor.
  static detran_utilities::SP<Quadrature>
    Create(size_t dim,
           size_t num_azimuths_octant,
           size_t num_space,
           size_t num_polar,
           std::string polar)
  {
    SP_quadrature
      p(new UniformCyclic(dim, num_azimuths_octant, num_space, num_polar, polar));
    return p;
  }

  ~UniformCyclic(){};

private:

};

} // end namespace detran_angle

#endif // detran_angle_UNIFORM_HH_

//---------------------------------------------------------------------------//
//              end of file UniformCyclic.hh
//---------------------------------------------------------------------------//
