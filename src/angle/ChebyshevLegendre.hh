//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ChebyshevLegendre.hh
 *  @brief ChebyshevLegendre class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_angle_CHEBYSHEVLEGENDRE_HH_
#define detran_angle_CHEBYSHEVLEGENDRE_HH_

#include "ProductQuadrature.hh"

namespace detran_angle
{

/**
 *  @class ChebyshevLegendre
 *  @brief Chebyshev-Legendre quadrature
 *
 *  The ChebyshevLegendre quadrature uses Gauss-Chebyshev quadrature
 *  over one azimuthal quadrant and Gauss-Legendre quadrature over
 *  the polar range.
 *
 *  Relevant database parameters:
 *    - quad_number_polar_octant    -- number polar angles per octant
 *    - quad_number_azimuth_octant  -- number azimuths per octant
 */
class ANGLE_EXPORT ChebyshevLegendre: public ProductQuadrature
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param dim    Dimension (2 or 3)
   *  @param na     Number of azimuths per octant
   *  @param np     Number of polar angles per octant
   */
  ChebyshevLegendre(const size_t dim, const size_t na, const size_t np);

  /// Virtual destructor
  virtual ~ChebyshevLegendre(){};

  /// SP constructor
  static SP_quadrature Create(const size_t dim,
                              const size_t na,
                              const size_t np);


};

} // end namespace detran_angle


#endif /* detran_angle_CHEBYSHEVLEGENDRE_HH_ */
