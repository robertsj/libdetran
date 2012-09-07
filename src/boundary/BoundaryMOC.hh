//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.hh
 * \brief  BoundaryMOC class definition
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYMOC_HH_
#define BOUNDARYMOC_HH_

#include "BoundaryBase.hh"
#include "BoundaryConditionMOC.hh"
#include "geometry/MeshMOC.hh"
#include "angle/QuadratureMOC.hh"

namespace detran
{

/*
 *  \class BoundaryMOC
 *  \brief Boundary flux container for MOC problems.
 *
 *  The method of characteristics solves the transport equation
 *  by sweeping along fixed tracks crossing the domain.  Tracks
 *  begin and end at a global boundary.  This class stores the
 *  angular flux for each track.
 *
 */

template <class D>
class BoundaryMOC : public BoundaryBase<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryBase<D>                     Base;
  typedef typename Base::SP_boundary          SP_boundary;
  typedef typename Base::SP_input             SP_input;
  typedef detran_geometry::Mesh               Mesh;
  typedef detran_geometry::MeshMOC            MeshMOC;
  typedef detran_utilities::SP<MeshMOC>       SP_mesh;
  typedef detran_angle::QuadratureMOC         QuadratureMOC;
  typedef detran_utilities::SP<QuadratureMOC> SP_quadrature;
  typedef BoundaryConditionMOC<D>             BC_T;
  typedef typename BC_T::SP_bc                SP_bc;
  typedef detran_utilities::size_t            size_t;
  typedef detran_utilities::vec_int           vec_int;
  typedef detran_utilities::vec2_int          vec2_int;
  typedef detran_utilities::vec3_int          vec3_int;
  typedef detran_utilities::vec2_dbl          vec2_dbl;
  typedef detran_utilities::vec3_dbl          vec3_dbl;
  typedef std::vector<vec3_dbl>               boundary_flux_type;

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
  BoundaryMOC(SP_input        input,
              SP_mesh         mesh,
              SP_quadrature   quadrature);

  /// SP Constructor.
  static SP_boundary
  Create(SP_input       input,
         SP_mesh        mesh,
         SP_quadrature  quadrature)
  {
    SP_boundary p(new BoundaryMOC(input, mesh, quadrature));
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
  void set_incident(const size_t g, double *v)
  {
    THROW("IMPLEMENT ME");
  }

  /*
   *  \brief Get the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void get_incident(const size_t g, double *v)
  {
    THROW("IMPLEMENT ME");
  }

  //-------------------------------------------------------------------------//
  // BOUNDARY FLUX ACCESS
  //-------------------------------------------------------------------------//

  /*
   *  \brief Const access to a boundary flux using cardinal indices.
   *
   *  This (and the mutable version) interface is for use
   *  in sweeping, where octants and angles are cycled.
   *
   *  \param    g     Energy group.
   *  \param    o     Octant index.
   *  \param    a     Angle index (within octant).
   *  \param    t     Track index for angle.
   *  \return         Constant reference to boundary flux.
   */
  const double&
  operator()(const size_t g, const size_t o, const size_t a,
             const size_t inout, const size_t t) const;

  /*
   *  \brief Mutable access to boundary flux using cardinal indices.
   */
  double&
  operator()(const size_t g, const size_t o, const size_t a,
             const size_t inout, const size_t t);

  /*
   *  \brief Assign the octant, angle, and track fed by an octant,
   *         angle, and track
   *
   *  To make the boundary condition as independent as
   *  possible, BoundaryMOC knows what to do with the
   *  cardinal angle within an octant (i.e. with polar implicit)
   *
   *  \param o1   Incident octant
   *  \param a1   Incident angle within octant
   *  \param t1   Incident track
   *  \param o2   Outgoing octant
   *  \param a2   Outgoing angle within octant
   *  \param t2   Outgoing track
   */
  void feed_into(const size_t o1, const size_t a1, const size_t t1,
                 size_t &o2, size_t &a2, size_t &t2);

  /*
   *  \brief Assign the octant, angle, and track feeding an octant,
   *         angle, and track
   *  \param o1   Outgoing octant
   *  \param a1   Outgoing angle within octant
   *  \param t1   Outgoing track
   *  \param o2   Incident octant
   *  \param a2   Incident angle within octant
   *  \param t2   Incident track
   */
  void feed_from(const size_t o1, const size_t a1, const size_t t1,
                 size_t &o2, size_t &a2, size_t &t2);

  /// Return vector of octant, azimuth, track triplets for a side.
  const vec2_int& side_indices(const size_t side) const
  {
    return d_side_index[side];
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

  /// MOC Quadrature
  SP_quadrature d_quadrature;

  /// Boundary flux [energy, angle, inout, track]
  boundary_flux_type d_boundary_flux;

  /// Vector of boundary conditions.
  std::vector<SP_bc> d_bc;

  /// d_feed_into[o0][a0][t0][o1 a1 t1] says track a0, t0 feeds into track a1, t1
  std::vector<vec3_int> d_feed_into;

  /// d_feed_from[a0][t0][a1 t1] says track a0, t0 feeds from track a1, t1
  std::vector<vec3_int> d_feed_from;

  /// d_side_index[side][o a t]
  vec3_int d_side_index;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Size the boundaries, etc.
  void initialize();

  /// Setup reflection indices.
  void setup_indices();

  /// Setup indices for an incident side.
  void setup_side_indices();

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "BoundaryMOC.i.hh"

#endif // BOUNDARYMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.hh
//---------------------------------------------------------------------------//
