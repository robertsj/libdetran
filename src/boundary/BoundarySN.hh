//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundarySN.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  BoundarySN class definition.
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYSN_HH_
#define BOUNDARYSN_HH_

#include "BoundaryBase.hh"
#include "BoundaryTraits.hh"
#include "BoundaryCondition.hh"
#include "angle/Quadrature.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class BoundarySN
 * \brief Boundary flux container for SN problems.
 *
 * Since arbitrary boundary functions are integral to the response matrix
 * method, it helps to have a really easy interface for handling boundary
 * information.  The Boundary classSN contains all boundary fluxes (incident
 * and outgoing) for each surface.  The exact type of a boundary flux
 * is templated. For example, a 1-D boundary flux (for an angle and group) is
 * just a double, whereas for a 3-D problem, it is a 2-D array.
 *
 * Keeping the fundamental type templated should allow this to be used in
 * all the SN stuff and potentially MOC applications. 
 *
 * All boundary fluxes are stored as follows: side->energy->angle->space^D
 *
 * Moreover, the angles are ordered in the same order as the quadrature,
 * i.e. octants->angles.
 *
 * To facilitate response generation somewhat, it is assumed that the
 * boundary fluxes are stored spatial in a way that allows fluxes to be
 * read from left to right, from bottom up on a face.  That means that
 * a sweeper can sweep over ii=0..I-1, jj= ..., selecte i, j, and k for
 * the octant (i.e. does it reverse?), but the boundaries are read using
 * the ii, jj, and kk.  What this allows for is that an incident boundary
 * condition can be inserted in a left to right, bottom up orientation
 * for any face.
 *
 * Unfortunately, the angular orientation isn't so easy.  Ideally, we'd
 * like phi=0..2-pi, theta=0..pi with respect to a surface, but that
 * would
 *
 */
//---------------------------------------------------------------------------//
template <class D>
class BoundarySN : public BoundaryBase<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryBase<D>                           Base;
  typedef typename Base::SP_boundary                SP_base;
  typedef detran_utilities::SP<BoundarySN<D> >      SP_boundary;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef typename BoundaryCondition<D>::SP_bc      SP_bc;
  typedef typename std::vector<boundary_flux_type>  vec_boundary_flux;
  typedef typename std::vector<vec_boundary_flux>   vec2_boundary_flux;
  typedef typename std::vector<vec2_boundary_flux>  vec3_boundary_flux;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::size_t                  size_t;

  using Base::IN;
  using Base::OUT;

  //-------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    quadrature  Angular quadrature.
   */
  BoundarySN(SP_input        input,
             SP_mesh         mesh,
             SP_quadrature   quadrature);

  /// SP Constructor
  static SP_base
  Create(SP_input         input,
         SP_mesh          mesh,
         SP_quadrature    quadrature)
  {
    SP_boundary p(new BoundarySN(input, mesh, quadrature));
    return p;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  \param  g   Group of current solve
   */
  void set(const size_t g);

  /*!
   *  \brief Update the boundaries for each sweep.
   *
   *  This updates all incident boundary fluxes
   *  using the current outgoing boundary fluxes
   *  in a group.  What happens in the update is
   *  a function of the underlying boundary
   *  condition.
   *
   *  \param  g   Group of current solve
   */
  void update(const size_t g);

  /*!
   *  \brief Update the boundaries for a single angle.
   *
   *  This is an alternative update that only updates
   *  the incident boundary flux for a particular
   *  angle.  When called within a sweep, this allows
   *  the most recent boundary fluxes to be used,
   *  in effect producing a Gauss-Seidel iteration.
   *
   *  \note This cannot be used for Krylov solvers.
   *
   *  \param  g   Group of current solve
   *  \param  o   Octant being swept
   *  \param  a   Angle within octant being swept
   */
  void update(const size_t g, const size_t o, const size_t a);

  /*!
   *  \brief Clear the boundary container for a group.
   *
   *  In some cases, a client might require homogeneous
   *  boundaries, perhaps after a fixed boundary has
   *  been used to construct a right hand side for a
   *  Krylov solve.
   *
   *  \param  g   Group of current solve
   */
  void clear(const size_t g);

  /*
   *  \brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is is support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void set_incident(const size_t g, double *v);

  /*
   *  \brief Get the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void get_incident(const size_t g, double *v);

  //-------------------------------------------------------------------------//
  // BOUNDARY FLUX ACCES
  //-------------------------------------------------------------------------//

  /*
   *  \brief Const access to a boundary flux using cardinal indices.
   *
   *  This (and the mutable version) interface is for use
   *  in sweeping, where octants and angles are cycled.
   *
   *  \param    side  Side index.
   *  \param    o     Octant index.
   *  \param    a     Angle index (within octant).
   *  \param    g     Energy group.
   *  \return         Constant reference to boundary flux.
   */
  const boundary_flux_type&
  operator()(const size_t side, const size_t o, const size_t a, const size_t g) const;

  /*
   *  \brief Mutable access to boundary flux using cardinal indices.
   */
  boundary_flux_type&
  operator()(const size_t side, const size_t o, const size_t a, const size_t g);

  /*
   *  \brief Const access to ordered incident flux.
   */
  const boundary_flux_type&
  incident(const size_t side, const size_t angle, const size_t g) const;

  /*
   *  \brief Mutable access to ordered incident boundary flux.
   */
  boundary_flux_type&
  incident(const size_t side, const size_t angle, const size_t g);

  /*
   *  \brief Const access to ordered outgoing flux.
   */
  const boundary_flux_type&
  outgoing(const size_t side, const size_t angle, const size_t g) const;

  /*
   *  \brief Mutable access to ordered outgoing boundary flux.
   */
  boundary_flux_type&
  outgoing(const size_t side, const size_t angle, const size_t g);

  /*!
   *  \brief  Map the local index to cardinal index.
   *
   *  In some cases, we need the boundary in its local order,
   *  e.g. left to right in space and angle.
   *
   *  \param  side    Side of interest
   *  \param  angle   Cardinal
   *  \param  inout
   *
   */
  size_t ordered_angle(const size_t side,
                       const size_t angle,
                       const size_t inout) const;

  //-------------------------------------------------------------------------//
  // GETTERS
  //-------------------------------------------------------------------------//

  /// Return the input.
  SP_input get_input() const
  {
    return d_input;
  }

  /// Return the mesh.
  SP_mesh get_mesh() const
  {
    return d_mesh;
  }

  /// Return the quadrature.
  SP_quadrature get_quadrature() const
  {
    return d_quadrature;
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Expose base class members.
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_number_groups;
  using Base::d_has_reflective;
  using Base::d_is_reflective;
  using Base::d_boundary_flux_size;

  /// Quadrature
  SP_quadrature d_quadrature;

  /// Boundary flux (side, energy, angle).(space^D)
  vec3_boundary_flux d_boundary_flux;

  /// Vector of boundary conditions.
  std::vector<SP_bc> d_bc;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Sizes the boundary flux containers.
  void initialize();

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "BoundarySN.i.hh"

#endif /* BOUNDARYSN_HH_ */

//---------------------------------------------------------------------------//
//              end of BoundarySN.hh
//---------------------------------------------------------------------------//
