//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   ChebyshevDPN.hh
 *  @brief  ChebyshevDPN
 *  @author Jeremy Roberts
 *  @date   Oct 11, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_CHEBYSHEVDPN_HH_
#define detran_angle_CHEBYSHEVDPN_HH_

#include "ProductQuadrature.hh"

namespace detran_angle
{

/**
 *  @class ChebyshevDPN
 *  @brief Chebyshev-Double PN quadrature
 *
 *  The ChebyshevDPN quadrature uses Gauss-Chebyshev quadrature
 *  over one azimuthal quadrant and DPN quadrature over
 *  the polar range.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant    -- number polar angles per octant
 *    - quad_number_azimuth_octant  -- number azimuths per octant
 */
class ANGLE_EXPORT ChebyshevDPN: public ProductQuadrature
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
  ChebyshevDPN(const size_t dim,
               const size_t na,
               const size_t np);

  /// Virtual destructor
  virtual ~ChebyshevDPN(){};

  /// SP constructor
  static SP_quadrature Create(const size_t dim,
                              const size_t na,
                              const size_t np);


};

} // end namespace detran_angle

#endif // detran_angle_CHEBYSHEVDPN_HH_

//---------------------------------------------------------------------------//
//              end of file ChebyshevDPN.hh
//---------------------------------------------------------------------------//
