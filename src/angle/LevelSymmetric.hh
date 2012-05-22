//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LevelSymmetric.hh
 * \author robertsj
 * \date   May 22, 2012
 * \brief  LevelSymmetric class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//
#ifndef LEVELSYMMETRIC_HH_
#define LEVELSYMMETRIC_HH_

// Detran
#include "Quadrature.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class LevelSymmetric
 * \brief 2D/3D Level-symmetric (LQn) quadrature class.
 *
 * Level symmetric quadratures are characterized by using the same set of
 * \f$ N/2 \f$ positive values of the direction cosines with respect to
 * each of the axes.  There are \f$ N(N+2)/8 \f$ ordinates per octant,
 * yielding \f$ N(N+2) \f$ directions for a 3D problem.
 *
 * Not all the direction cosines are independent; for a given \f$ \mu_n \f$,
 * one has a single degree of freedom, e.g. one can further choose only
 * \f$ \eta_n \f$.  The final cosine is of course defined such that
 * the sum in quadrature is unity.
 *
 * An unfortunate aspect of LQn is that negative weights appear for
 * \f$ N = 20 \f$.  As an alternative, see the uniform and equal
 * weight quadrature set (UEn) \ref Uniform_Equal_3D, which yields positive
 * weights for arbitrarily high \f$ N \f$.
 *
 * \param order   quadrature order
 * \param N       spatial dimension (default 3)
 *
 */
//===========================================================================//
class LevelSymmetric : public Quadrature
{

public:

  typedef SP<LevelSymmetric>  SP_quadrature;
  typedef SP<Quadrature>      SP_base;

  /*!
   *  \brief Constructor.
   *
   *  \param    order       Quadrature order.
   *  \param    dim         Problem dimension
   */
  LevelSymmetric(int order, int dim);

  /*!
   *  \brief SP Constructor.
   */
  static SP<Quadrature> Create(int order, int dim)
  {
    SP_quadrature p;
    p = new LevelSymmetric(order, dim);
    return p;
  }

  void display() const;

private:

  /// Pre-computed tables of abscissa and weights.
  void set_quad_values(vec_dbl &att,
                       vec_dbl &wtt);

};

} // end namespace detran

#endif /* LEVELSYMMETRIC_HH_ */

//---------------------------------------------------------------------------//
//              end of Level_Symmetric_3D.hh
//---------------------------------------------------------------------------//
