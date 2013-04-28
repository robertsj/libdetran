//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  UniformEqual.hh
 *  @brief UniformEqual class definition
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_UNIFORMEQUAL_HH_
#define detran_angle_UNIFORMEQUAL_HH_

#include "Quadrature.hh"

namespace detran_angle
{

//---------------------------------------------------------------------------//
/**
 *  @class UniformEqual
 *  @brief 2D/3D Uniform, Equal Weight (UEn) quadrature class.
 *
 * As mentioned in \ref LevelSymmetric, a fundamental problem inherent
 * to LQn quadrature is presence of negative weights at high order.  These
 * weights produce unphysical solutions (and may inhibit convergence).  For
 * problems where an increased quadrature order (i.e. more angles) is
 * required to study convergence or simply to get better answers, we
 * require an arbitrarily high order, positive weight quadrature.
 *
 * Here, we implement the "uniform, equal weight" (UEn) quadrature of Carew
 * and Zamonsky.  The basic idea is to choose uniform azimuthal divisions
 * and uniform polar cosines.  The result is a product quadrature that
 * can be extended on-the-fly to arbitrary numbers of angles.
 *
 * Here, the quadrature order defines the number of polar angles.
 * We take the number of azimuthal angles to be twice this
 * number.  Hence, the total number of angles is
 * twice the number of polar angles squared.
 *
 * Note, the angles are stored with polar as the inner index.  This is
 * requested for all product quadratures.
 *
 * \refs
 * - Carew, J. and Zamonsky. G., <em>Nuclear Science and Engineering</em>
 *     <b>131</b>, 199-207 (1999).
 *
 */
//---------------------------------------------------------------------------//

class ANGLE_EXPORT UniformEqual : public Quadrature
{

public:

  //---------------------------------------------------------------------------//
  // TYPEDEFS
  //---------------------------------------------------------------------------//

  typedef detran_utilities::SP<UniformEqual>  SP_quadrature;
  typedef detran_utilities::SP<Quadrature>    SP_base;

  //---------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //---------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    order       Quadrature order.
   *  @param    dim         Problem dimension
   */
  UniformEqual(size_t order, size_t dim);

  /// SP constructor
  static detran_utilities::SP<Quadrature>
  Create(size_t order, size_t dim)
  {
    SP_quadrature p(new UniformEqual(order, dim));
    return p;
  }

};

} // end namespace detran_angle

#endif /* detran_angle_UNIFORMEQUAL_HH_ */

//---------------------------------------------------------------------------//
//              end of UniformEqual.hh
//---------------------------------------------------------------------------//
