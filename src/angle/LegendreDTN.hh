//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LegendreDTN.hh
 *  @brief  LegendreDTN
 *  @author Jeremy Roberts
 *  @date   Jan 23, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_LEGENDREDTN_HH_
#define detran_angle_LEGENDREDTN_HH_

#include "ProductQuadrature.hh"

namespace detran_angle
{

/**
 *  @class LegendreDTN
 *  @brief Legendre-Double TN quadrature
 *
 *  The LegendreDTN quadrature uses Gauss-Legendre quadrature
 *  over one azimuthal quadrant and DTN quadrature over
 *  the polar range.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant    -- number polar angles per octant
 *    - quad_number_azimuth_octant  -- number azimuths per octant
 */
class LegendreDTN: public ProductQuadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param dim    Dimension (2 or 3)
   *  @param na     Number of azimuths per octant
   *  @param np     Number of polar angles per octant
   *  @param name   Quadrature name
   */
  LegendreDTN(const size_t dim,
               const size_t na,
               const size_t np);

  /// Virtual destructor
  virtual ~LegendreDTN(){};

  /// SP constructor
  static SP_quadrature Create(const size_t dim,
                              const size_t na,
                              const size_t np);


};

} // end namespace detran_angle

#endif // detran_angle_LEGENDREDTN_HH_

//---------------------------------------------------------------------------//
//              end of file LegendreDTN.hh
//---------------------------------------------------------------------------//
