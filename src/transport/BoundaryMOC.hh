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

// Detran
#include "BoundaryBase.hh"
#include "BoundaryConditionMOC.hh"
#include "MeshMOC.hh"
#include "QuadratureMOC.hh"

// Utilities
#include "DBC.hh"

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

  typedef SP<BoundaryMOC>               SP_boundary;
  typedef InputDB::SP_input             SP_input;
  typedef SP<MeshMOC>                   SP_mesh;
  typedef SP<QuadratureMOC>             SP_quadrature;
  typedef BoundaryConditionMOC<D>       BC_T;
  typedef typename BC_T::SP_bc          SP_bc;
  typedef BoundaryBase<D>               Base;
  typedef typename Base::SP_boundary    SP_base;
  typedef std::vector<vec3_dbl>         boundary_flux_type;

  using Base::IN;
  using Base::OUT;

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
  static SP<detran::BoundaryBase<D> >
  Create(SP<detran::InputDB>       input,
         SP<detran::Mesh>          mesh,
         SP<detran::Quadrature>    quadrature)
  {
    SP_boundary p(new BoundaryMOC(input, mesh, quadrature));
    return p;
  }

  /// \name Inherited Interface
  /// \{

  /*!
   *  \brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  \param  g   Group of current solve
   */
  void set(int g);

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
  void update(int g);

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
  void update(int g, int o, int a);

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
  void clear(int g);

  /*
   *  \brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is is support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void set_incident(int g, double *v)
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
  void get_incident(int g, double *v)
  {
    THROW("IMPLEMENT ME");
  }

  /// \}

  /// \name Track Boundary Flux Access
  /// \{

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
  operator()(u_int g, u_int o, u_int a, u_int inout, u_int t) const;

  /*
   *  \brief Mutable access to boundary flux using cardinal indices.
   */
  double&
  operator()(u_int g, u_int o, u_int a, u_int inout, u_int t);

  /// \}

  /// \name Indexing
  /// \{

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
  void feed_into(const u_int o1, const u_int a1, const u_int t1,
                 u_int &o2, u_int &a2, u_int &t2);

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
  void feed_from(const u_int o1, const u_int a1, const u_int t1,
                 u_int &o2, u_int &a2, u_int &t2);

  /// Return vector of octant, azimuth, track triplets for a side.
  const vec2_int& side_indices(int side)
  {
    return d_side_index[side];
  }

  /// \}

  /// DBC function
  bool is_valid() const
  {
    return true;
  }

private:

  // Expose base class members.
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_number_groups;
  using Base::d_has_reflective;
  using Base::d_is_reflective;
  using Base::d_boundary_flux_size;

  /// \name Private Data
  /// \{

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

  /// \}

  /// \name Implementation
  /// \{

  /// Size the boundaries, etc.
  void initialize();

  /// Setup reflection indices.
  void setup_indices();

  /// Setup indices for an incident side.
  void setup_side_indices();

  /// \}

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
