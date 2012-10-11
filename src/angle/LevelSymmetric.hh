//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   LevelSymmetric.hh
 *  @author robertsj
 *  @date   May 22, 2012
 *  @brief  LevelSymmetric class definition.
 */
//---------------------------------------------------------------------------//
#ifndef LEVELSYMMETRIC_HH_
#define LEVELSYMMETRIC_HH_

#include "Quadrature.hh"

namespace detran_angle
{

/**
 *  @class LevelSymmetric
 *  @brief 2D/3D Level-symmetric (LQn) quadrature class.
 *
 *  Level symmetric quadratures are characterized by using the same set of
 *  \f$ N/2 \f$ positive values of the direction cosines with respect to
 *  each of the axes.  There are \f$ N(N+2)/8 \f$ ordinates per octant,
 *  yielding \f$ N(N+2) \f$ directions for a 3D problem.
 *
 *  Not all the direction cosines are independent; for a given \f$ \mu_n \f$,
 *  one has a single degree of freedom, e.g. one can further choose only
 *  \f$ \eta_n \f$.  The final cosine is of course defined such that
 *  the sum in quadrature is unity.
 *
 *  An unfortunate aspect of LQn is that negative weights appear for
 *  \f$ N = 20 \f$.  As an alternative, see the uniform and equal
 *  weight quadrature set (UEn) \ref UniformEqual, which yields positive
 *  weights for arbitrarily high \f$ N \f$.
 *
 */
class LevelSymmetric : public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  \brief Constructor.
   *  \param    order       Quadrature order.
   *  \param    dim         Problem dimension
   */
  LevelSymmetric(size_t order, size_t dim);

  /**
   *  \brief SP Constructor.
   */
  static SP_quadrature Create(size_t order, size_t dim)
  {
    SP_quadrature p;
    p = new LevelSymmetric(order, dim);
    return p;
  }

private:

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Pre-computed tables of abscissa and weights.
  void set_quad_values(const size_t order,
                       vec_dbl &att,
                       vec_dbl &wtt);

};

} // end namespace detran_angle

#endif /* LEVELSYMMETRIC_HH_ */

//---------------------------------------------------------------------------//
//              end of Level_Symmetric_3D.hh
//---------------------------------------------------------------------------//
