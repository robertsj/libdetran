//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   MomentToDiscrete.hh
 * \author Jeremy Roberts
 * \date   Jul 1, 2011
 * \brief  MomentToDiscrete class definition.
 * \note   Copyright (C) 2011 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef MOMENT_TO_DISCRETE_HH_
#define MOMENT_TO_DISCRETE_HH_

// Detran
#include "Quadrature.hh"

// Utilities
#include "Definitions.hh"
#include "SP.hh"
#include "Traits.hh"

// System
#include <vector>

namespace detran
{

//===========================================================================//
/*!
 * \class MomentToDiscrete
 * \brief The discrete ordinates moment-to-discrete operator.
 *
 * This class defines the operator
 * \f[
   \mathbf{M} =
   \left(\begin{array}{llllllll}
     \frac{1}{4\pi}Y^{0}_{0}(\Omega_1)         &
     \frac{3}{4\pi}Y^{-1}_{1}(\Omega_1)        &
     \frac{3}{4\pi}Y^{0}_{1}(\Omega_1)         &
     \frac{3}{4\pi}Y^{1}_{1}(\Omega_1)         &
     \frac{5}{4\pi}Y^{-2}_{2}(\Omega_1)        &
     \ldots                                     &
     \frac{2L+1}{4\pi}Y^{L-1}_{L}(\Omega_1)    &
     \frac{2L+1}{4\pi}Y^{L}_{L}(\Omega_1)      \\
     \frac{1}{4\pi}Y^{0}_{0}(\Omega_2)         &
     \frac{3}{4\pi}Y^{-1}_{1}(\Omega_2)        &
     \frac{3}{4\pi}Y^{0}_{1}(\Omega_1)         &
     \frac{3}{4\pi}Y^{1}_{1}(\Omega_2)         &
     \frac{5}{4\pi}Y^{-2}_{2}(\Omega_2)        &
     \ldots                                     &
     \frac{2L+1}{4\pi}Y^{L-1}_{L}(\Omega_2)    &
     \frac{2L+1}{4\pi}Y^{L}_{L}(\Omega_2)      \\
     \vdots & \vdots & \vdots & \vdots &
     \vdots & & \vdots & \vdots                 \\
     \frac{1}{4\pi}Y^{0}_{0}(\Omega_{N_n})     &
     \frac{3}{4\pi}Y^{-1}_{1}(\Omega_{N_n})    &
     \frac{3}{4\pi}Y^{0}_{1}(\Omega_{N_n})     &
     \frac{3}{4\pi}Y^{1}_{1}(\Omega_{N_n})     &
     \frac{5}{4\pi}Y^{-2}_{2}(\Omega_{N_n})    &
     \ldots                                     &
     \frac{2L+1}{4\pi}Y^{L-1}_{L}(\Omega_{N_n})&
     \frac{2L+1}{4\pi}Y^{L}_{L}(\Omega_{N_n})  \\
   \end{array}\right) \, .
 * \f]
 *
 * Then, for instance, the angular flux in a particular direction
 * for a particular cell and group can be approximated by the dot
 * product of the corresponding row of \f$\mathbf{M}\f$ with the
 * vector of flux moments for that cell and group, i.e.
 * \f[
 *    \psi_{i,n} \approx \frac{1}{4\pi} Y^{0}_{0}(\Omega_n)\phi^{0}_{0}    +
 *                       \frac{3}{4\pi} Y^{-1}_{1}(\Omega_n)\phi^{-1}_{1}  +
 *                       \ldots                                            +
 *                       \frac{2L+1}{4\pi} Y^{L}_{L}(\Omega_n)\phi^L_L  \, ,
 * \f]
 * for cell \f$i\f$ and angle \f$n\f$.  For a purely isotropic flux, we
 * have \f$\phi=\phi^{0}_{0}\f$ and \f$\psi_n = \phi/4\pi\f$ as we expect
 * (since \f$ Y^{0}_{0} = 1 \f$).
 *
 * For anisotropic scattering  of order \f$L\f$ in 3-d problems,
 * the number of columns in \f$\mathbf{M}\f$ is \f$ N_L =(L+1)^2\f$.  In
 * 2-d problems, there is no variation in the \f$z\f$ direction, and hence
 *
 * For 1-d problems and scattering order \f$L\f$, this reduces to
 * \f$ N_L = (L+1) \f$ and a subsequent simplification of \f$\mathbf{M}\f$
 * into terms only of Legendre polynomials (and a normalization of
 * \f$ (2l+1)/2 \f$ instead of \f$(2l+1)/4\pi \f$.
 *
 * The number of rows \f$ N_n \f$ is determined by the number of angles,
 * which is defined by the quadrature and its order.
 *
 * The moments of an expanded function are organized as a row of
 * \f$\mathbf{M}\f$, in an array, and one can easily index into
 * that array for a given \f$l\f$ and \f$m\f$ using the
 * \ref Moments::index of appropriated dimension.
 *
 * \tparam D dimension
 *
 * \sa Spherical_Harmonics, Moments
 *
 */
/*!
 * \example test/testMomentToDiscrete.cc
 *
 * Test of MomentToDiscrete.
 */
//===========================================================================//

template <class D>
class MomentToDiscrete
{

public:

  /// \name Useful Typedefs
  //\{
  typedef SP<MomentToDiscrete>              SP_MtoD;
  typedef unsigned int                      size_type;
  typedef vec_dbl                           M_Row;
  typedef std::vector<M_Row>                Operator_M;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  //\}

  /*!
   * \brief Constructor.
   *
   * \param legendre_order   Legendre order of moments expansion.
   *
   * Note, this will be greater than or equal to the flux moments
   * order, since we can optionally keep higher order moments.
   */
  explicit MomentToDiscrete(const size_type legendre_order);

  /// SP contructor
  static SP_MtoD Create(const size_type o, SP_quadrature q)
  {
    SP_MtoD m(new MomentToDiscrete(o));
    m->build(q);
    return m;
  }

  /*!
   * \brief Build the moments-to-discrete operator.
   *
   * Keeping the actual construction outside the constructor allows
   * us to rebuild the operator for different angular solves using
   * the same spatial grid. This is useful for coupled forward and
   * adjoint solves, compact testing of quadrature sets, and potential
   * angular multigrid schemes.
   *
   * \param     angularmesh     Angular mesh smart pointer.
   */
  void build(SP_quadrature quadrature);




  /// \name Accessors
  //\{

  /*!
   * \brief Return an element from \f$ M\f$.
   *
   * \param     angle       Angle i.e. row.
   * \param     moment      Moment cardinal index, i.e. column.
   * \return                Element of operator.
   */
  const double& operator()(const int angle, const int moment) const;

  /*!
   * \brief Return an element from \f$ M\f$.
   *
   * \param     angle       Angle i.e. row.
   * \param     l           Legendre degree.
   * \param     m           Legendre order.
   * \return                Element of operator.
   */
  const double& operator()(const int angle, const int l, const int m) const;


  /*!
   * \brief Return an element from \f$ M\f$.
   *
   * \param     o           Octant index.
   * \param     a           Angle index (within octant).
   * \param     l           Legendre degree.
   * \param     m           Legendre order.
   * \return                Element of operator.
   */
  const double& operator()(const int o, const int a,
                           const int l, const int m) const;

  /*!
   * \brief Return a row of the operator.
   *
   * \param     angle       Angle i.e. row index.
   * \return                A row.
   */
  const M_Row& get_row(const int angle) const
  {
    Require(angle > 0);
    Require(angle < d_number_angles);
    return d_M[angle];
  }

  /// Return number of moments (length of row in \f$\mathbf{M}\f$).
  size_type row_size() const
  {
    return d_number_moments;
  }

  /// Return number of angles (length of column in \f$\mathbf{M}\f$).
  size_type column_size() const
  {
    return d_number_angles;
  }

  /// Return Legendre expansion order.
  size_type legendre_order() const
  {
    return d_legendre_order;
  }


  // \}

private:

  // >>> DATA

  /// Legendre order of anisotropic scattering.
  const size_type d_legendre_order;

  /// Number of angular moments.
  const size_type d_number_moments;

  /// Angular mesh.
  SP_quadrature d_quadrature;

  /// Number of angles (not const so that we can change angular meshes)
  size_type d_number_angles;

  /// Moments-to-discrete operator \f$\mathbf{M}\f$.
  Operator_M d_M;

  /// \name Implementation
  //\{

  /*!
   * \brief Calculate one row of \f$\mathbf{M}\f$.
   *
   * \param o   Octant index.
   * \param a   Angle index.
   */
  void calc_row(const int o, const int a);

  //\}

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "MomentToDiscrete.i.hh"

#endif /* MOMENT_TO_DISCRETE_HH_ */

//---------------------------------------------------------------------------//
//              end of MomentToDiscrete.hh
//---------------------------------------------------------------------------//
