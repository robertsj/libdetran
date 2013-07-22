//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  ContinuousOrthogonalBasis.hh
 *  @brief ContinuousOrthogonalBasis class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_orthog_CONTINUOUSORTHOGONALBASIS_HH_
#define detran_orthog_CONTINUOUSORTHOGONALBASIS_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class ContinuousOrthogonalBasis
 *  @brief Represents a continuous basis on an arbitrary domain.
 *
 *  See @ref OrthogonalBasis for some additional background.
 *
 *  For continuous functions defined on @f$ x \in [a, b] @f$,
 *  an expansion in an orthogonal basis defined on the same
 *  domain is defined
 *  @f[
 *      f =  \sum^{\infty}_{l=0} \int^a_b P_l(x) a_l^{-1} \tilde{f}_l \, ,
 *  @f]
 *  where the expansion coefficients (i.e. the transformed function) is
 *  defined
 *  @f[
 *      \tilde{f}_l =  \int^a_b P_l(x) w_l(x) f(x) dx \, .
 *  @f]
 *  On a discretized domain, we know @f$ f_i = f(x_i) @f$ at discrete
 *  points @f$ x_i @f$.  (In finite element analysis, we
 *  would actually have a continuous representation of functions,
 *  but Detran uses no finite element discretizations yet).
 *
 *  Hence, the integrals must be approximated.  The easiest
 *  way to do this for arbitrary meshes is via the midpoint
 *  rule, which while second order is actually a perfect fit
 *  for the discretizations used in Detran as they all assume
 *  a constant value over the mesh.  Then
 *  @f[
 *      \tilde{f}_l \approx
 *        \sum_{i=0}^{N} P_l(x_i) w_l(x_i) f(x_i) \Delta_i \, .
 *  @f]
 *  In general, any quadrature may be used.
 */
class ORTHOG_EXPORT ContinuousOrthogonalBasis: public OrthogonalBasis
{

public:

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// @brief Constructor
  ContinuousOrthogonalBasis(const Parameters &p);

  /// Pure virtual destructor
  virtual ~ContinuousOrthogonalBasis() = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Grid on which basis is defined
  vec_dbl d_x;
  /// Quadrature weights
  vec_dbl d_qw;
  /// Lower bound
  double d_lower_bound;
  /// Upper bound
  double d_upper_bound;

};

} // end namespace detran_orthog

#endif // detran_orthog_CONTINUOUSORTHOGONALBASIS_HH_

//----------------------------------------------------------------------------//
//              end of file ContinuousOrthogonalBasis.hh
//----------------------------------------------------------------------------//
