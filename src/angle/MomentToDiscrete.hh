//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MomentToDiscrete.hh
 *  @author Jeremy Roberts
 *  @date   Jul 1, 2011
 *  @brief  MomentToDiscrete class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_MOMENT_TO_DISCRETE_HH_
#define detran_angle_MOMENT_TO_DISCRETE_HH_

#include "Quadrature.hh"
#include "MomentIndexer.hh"
#include "utilities/Definitions.hh"
#include "utilities/SP.hh"
#include <vector>

namespace detran_angle
{

//---------------------------------------------------------------------------//
/**
 *  @class MomentToDiscrete
 *  @brief Converts moment-valued unknowns to discrete angle values
 *
 * This class defines the operator
 *  @f[
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
 *  @f]
 *
 * Then, for instance, the angular flux in a particular direction
 * for a particular cell and group can be approximated by the dot
 * product of the corresponding row of \f$\mathbf{M}\f$ with the
 * vector of flux moments for that cell and group, i.e.
 *  @f[
 *    \psi_{i,n} \approx \frac{1}{4\pi} Y^{0}_{0}(\Omega_n)\phi^{0}_{0}    +
 *                       \frac{3}{4\pi} Y^{-1}_{1}(\Omega_n)\phi^{-1}_{1}  +
 *                       \ldots                                            +
 *                       \frac{2L+1}{4\pi} Y^{L}_{L}(\Omega_n)\phi^L_L  \, ,
 *  @f]
 * for cell \f$i\f$ and angle \f$n\f$.  For a purely isotropic flux, we
 * have \f$\phi=\phi^{0}_{0}\f$ and \f$\psi_n = \phi/4\pi\f$ as we expect
 * (since \f$ Y^{0}_{0} = 1 \f$).
 *
 * For anisotropic scattering  of order \f$L\f$ in 3-d problems,
 * the number of columns in \f$\mathbf{M}\f$ is \f$ N_L =(L+1)^2\f$.  In
 * 2-d problems, there is no variation in the \f$z\f$ direction, and hence
 * all values with odd polar order vanish.
 *
 * For 1-d problems and scattering order \f$L\f$, this reduces to
 *  @f$ N_L = (L+1) \f$ and a subsequent simplification of \f$\mathbf{M}\f$
 * into terms only of Legendre polynomials (and a normalization of
 *  @f$ (2l+1)/2 \f$ instead of \f$(2l+1)/4\pi \f$.
 *
 * The number of rows \f$ N_n \f$ is determined by the number of angles,
 * which is defined by the quadrature and its order.
 *
 * The moments of an expanded function are organized as a row of
 *  @f$\mathbf{M}\f$, in an array, and one can easily index into
 * that array for a given \f$l\f$ and \f$m\f$ using the
 *  @ref Moments::index of appropriated dimension.
 *
 *  @sa Spherical_Harmonics
 *
 */
/**
 *  @example test/test_MomentToDiscrete.cc
 *
 * Test of MomentToDiscrete.
 */
//---------------------------------------------------------------------------//
class ANGLE_EXPORT MomentToDiscrete
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<MomentToDiscrete>  SP_MtoD;
  typedef MomentIndexer::SP_momentindexer         SP_momentindexer;
  typedef Quadrature::SP_quadrature               SP_quadrature;
  typedef detran_utilities::size_t                size_t;
  typedef detran_utilities::vec_dbl               M_Row;
  typedef std::vector<M_Row>                      Operator_M;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param indexer   Indexer for spherical harmonic orders
   */
  explicit MomentToDiscrete(SP_momentindexer indexer);

  /// SP contructor
  static SP_MtoD Create(SP_momentindexer indexer, SP_quadrature q)
  {
    SP_MtoD m(new MomentToDiscrete(indexer));
    m->build(q);
    return m;
  }

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /**
   *  \brief Build the moments-to-discrete operator.
   *
   *  Keeping the actual construction outside the constructor allows
   *  us to rebuild the operator for different angular solves using
   *  the same spatial grid. This is useful for coupled forward and
   *  adjoint solves, compact testing of quadrature sets, and potential
   *  angular multigrid schemes.
   *
   *  \param     q     Pointer to quadrature
   */
  void build(SP_quadrature q);

  /**
   *  @brief Return an element from \f$ M\f$.
   *
   *  @param     angle       Angle i.e. row.
   *  @param     moment      Moment cardinal index, i.e. column.
   *  @return                Element of operator.
   */
  const double& operator()(const size_t angle, const size_t moment) const;

  /**
   *  @brief Return an element from \f$ M\f$.
   *
   *  @param     angle       Angle i.e. row.
   *  @param     l           Legendre degree.
   *  @param     m           Legendre order.
   *  @return                Element of operator.
   */
  const double& operator()(const size_t angle,
                           const size_t l,
                           const size_t m) const;

  /**
   *  @brief Return an element from \f$ M\f$.
   *
   *  @param     o           Octant index.
   *  @param     a           Angle index (within octant).
   *  @param     l           Legendre degree.
   *  @param     m           Legendre order.
   *  @return                Element of operator.
   */
  const double& operator()(const size_t o, const size_t a,
                           const size_t l, const size_t m) const;

  /**
   *  @brief Return a row of the operator.
   *
   *  @param     angle       Angle i.e. row index.
   *  @return                A row.
   */
  const M_Row& get_row(const size_t angle) const
  {
    Require(angle < d_number_angles);
    return d_M[angle];
  }

  /// Return number of moments (length of row in \f$\mathbf{M}\f$).
  size_t row_size() const
  {
    return d_number_moments;
  }

  /// Return number of angles (length of column in \f$\mathbf{M}\f$).
  size_t column_size() const
  {
    return d_number_angles;
  }

  /// Return Legendre expansion order.
  size_t legendre_order() const
  {
    return d_legendre_order;
  }


private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Moment indexer
  SP_momentindexer d_indexer;
  /// Legendre order of anisotropic scattering.
  size_t d_legendre_order;
  /// Number of angular moments.
  size_t d_number_moments;
  /// Angular mesh.
  SP_quadrature d_quadrature;
  /// Number of angles (not const so that we can change angular meshes)
  size_t d_number_angles;
  /// Moments-to-discrete operator \f$\mathbf{M}\f$.
  Operator_M d_M;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /**
   *  @brief Calculate one row of \f$\mathbf{M}\f$.
   *
   *  @param o   Octant index.
   *  @param a   Angle index.
   */
  void calc_row_1d(const size_t o, const size_t a);
  void calc_row_2d(const size_t o, const size_t a);
  void calc_row_3d(const size_t o, const size_t a);

};

ANGLE_TEMPLATE_EXPORT(detran_utilities::SP<MomentToDiscrete>)

} // end namespace detran_angle

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "MomentToDiscrete.i.hh"

#endif /* detran_angle_MOMENT_TO_DISCRETE_HH_ */

//---------------------------------------------------------------------------//
//              end of MomentToDiscrete.hh
//---------------------------------------------------------------------------//
