//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  CurrentTally.hh
 *  @brief CurrentTally class definition
 *  @note  Copyright (C) Jeremy Roberts 2013
 */
//----------------------------------------------------------------------------//

#ifndef detran_CURRENTTALLY_HH_
#define detran_CURRENTTALLY_HH_

#include "transport/BoundaryTally.hh"

namespace detran
{

/**
 *  @class CurrentTally
 *  @brief Records partial currents through coarse mesh surfaces.
 *
 *  To illustrate, consider the leakage term of the transport
 *  equation in one coarse cell:
 *
 *  \f[
 *      \Big ( \mu  \frac{\partial}{\partial x}
 *           + \eta \frac{\partial}{\partial y}
 *           + \xi  \frac{\partial}{\partial z}
 *      \Big ) \psi(x, \mu, \eta, \xi) \, .
 *  \f]
 *
 *  Assuming a cell of volume \f$ V = \Delta_x \Delta_y \Delta_z \f$, we
 *  integrate the first term over the cell to get
 *
 *  \f[
 *  \begin{split}
 *      & \mu \int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz
 *         \Bigg ( \frac{\partial \psi}{\partial x} \Big |_{x=\Delta_x}
 *                -\frac{\partial \psi}{\partial x} \Big |_{x=0} \Bigg ) =
 *      \int^{\Delta_y}_{0} \int^{\Delta_z}_{0} dy dz
 *        \Big ( J_x(\Delta_x, y, z) - J_x(0, y, z) \Big ) \, ,
 *  \end{split}
 *  \f]
 *
 *  where \f$ J_x \f$ is the \f$x\f$-directed partial current.  Similar
 *  terms can be found for the other directions.
 *
 *  In 2D and 3D cases, the current must be integrated in space, which
 *  is done using a simple mid-point rule consistent with the underlying
 *  spatial discretizations currently available in Detran.  If higher
 *  order methods are implemented, equation-dependent currents would
 *  be required.
 *
 */
/**
 *  @example transport/test/test_CurrentTally.cc
 *
 *  Test of CurrentTally.
 */

template <class D>
class CurrentTally: public BoundaryTally<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryTally<D>                              Base;
  typedef detran_utilities::SP<CurrentTally>            SP_currenttally;
  typedef CoarseMesh::SP_coarsemesh                     SP_coarsemesh;
  typedef detran_angle::Quadrature::SP_quadrature       SP_quadrature;
  typedef CoarseMesh::SP_mesh                           SP_mesh;
  typedef typename EquationTraits<D>::face_flux_type    face_flux_type;
  typedef detran_utilities::size_t                      size_t;
  typedef detran_utilities::vec_int                     vec_int;
  typedef detran_utilities::vec2_int                    vec2_int;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::vec2_dbl                    vec2_dbl;
  typedef detran_utilities::vec3_dbl                    vec3_dbl;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mesh
   *  @param quadrature
   *  @param number_groups
   */
  CurrentTally(SP_coarsemesh mesh,
               SP_quadrature quadrature,
               const size_t  number_groups);

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE --- ALL TALLIES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /**
   *  @brief Add angular flux to the current tally
   *
   *  @param  i           x mesh index
   *  @param  j           y mesh index
   *  @param  k           z mesh index
   *  @param  g           group index
   *  @param  o           octant
   *  @param  a           angle within octant
   *  @param  psi         edge angular flux
   */
  void tally(const size_t i,
             const size_t j,
             const size_t k,
             const size_t g,
             const size_t o,
             const size_t a,
             const face_flux_type psi);

  /**
   *  @brief Add angular flux to the tally for a single incident direction
   *
   *  This is used to sweep over incident boundary containers so that
   *  the incident conditions can be avoided for the sweep.
   *
   *  @param  i           x mesh index
   *  @param  j           y mesh index
   *  @param  k           z mesh index
   *  @param  g           group index
   *  @param  o           octant
   *  @param  a           angle within octant
   *  @param  d           axis index for the incident flux
   *  @param  psi         edge angular flux
   */
  void tally(const size_t i,
             const size_t j,
             const size_t k,
             const size_t g,
             const size_t o,
             const size_t a,
             const size_t d,
             const double psi);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /**
   *  @brief Get the partial current from a surface and sense
   *  @param  i       x mesh index
   *  @param  j       y mesh index
   *  @param  k       z mesh index
   *  @param  g       group index
   *  @param  axis    0, 1, or 2 for x, y, or z
   *  @param  sense   true for positive (e.g. +x)
   */
  double partial_current(const size_t i,
                         const size_t j,
                         const size_t k,
                         const size_t g,
                         const size_t axis,
                         const size_t sense);

  /// Print all the partial currents (for debugging)
  void display();

  /// Reset a group
  void reset(const size_t group);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  // Expose base class members
  using Base::d_coarsemesh;
  using Base::d_quadrature;
  using Base::d_number_groups;
  using Base::d_coarse_edge_flag;
  using Base::d_octant_shift;
  using Base::d_perpendicular_index;
  using Base::index;

  /// Partial currents [dimension][group][sense][index]
  std::vector<vec3_dbl> d_partial_current;

};

} // end namespace detran

//----------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//----------------------------------------------------------------------------//

#include "CurrentTally.i.hh"

#endif // detran_CURRENTTALLY_HH_

//----------------------------------------------------------------------------//
//              end of file CurrentTally.hh
//----------------------------------------------------------------------------//
