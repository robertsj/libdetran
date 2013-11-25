//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryMOC.hh
 *  @brief  BoundaryMOC class definition
 *  @author Jeremy Roberts
 *  @date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYMOC_HH_
#define detran_BOUNDARYMOC_HH_

#include "boundary/BoundaryBase.hh"
#include "boundary/BoundaryConditionMOC.hh"
#include "geometry/Mesh.hh"
#include "angle/ProductQuadrature.hh"

namespace detran
{

/**
 *  @class BoundaryMOC
 *  @brief Boundary flux container for MOC problems.
 *
 *  The method of characteristics solves the transport equation
 *  by sweeping along fixed tracks crossing the domain.  Tracks
 *  begin and end at a global boundary.  This class stores the
 *  angular flux for each track.
 */

template <class D>
class BoundaryMOC : public BoundaryBase<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryBase<D>                       Base;
  typedef typename Base::SP_boundary            SP_base;
  typedef detran_utilities::SP<BoundaryMOC<D> > SP_boundary;
  typedef typename Base::SP_input               SP_input;
  typedef detran_geometry::Mesh                 Mesh;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef detran_angle::ProductQuadrature       QuadratureMOC;
  typedef detran_utilities::SP<QuadratureMOC>   SP_quadrature;
  typedef BoundaryConditionMOC<D>               BC_T;
  typedef typename BC_T::SP_bc                  SP_bc;
  typedef detran_utilities::size_t              size_t;
  typedef detran_utilities::vec_int             vec_int;
  typedef detran_utilities::vec2_int            vec2_int;
  typedef detran_utilities::vec3_int            vec3_int;
  typedef detran_utilities::vec2_dbl            vec2_dbl;
  typedef detran_utilities::vec3_dbl            vec3_dbl;
  typedef std::vector<vec3_dbl>                 bf_type;
  typedef D                                     D_T;

  using Base::IN;
  using Base::OUT;

  //-------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    input       User input database.
   *  @param    mesh        Cartesian mesh.
   *  @param    quadrature  Angular quadrature.
   */
  BoundaryMOC(SP_input        input,
              SP_mesh         mesh,
              SP_quadrature   quadrature);

  /// SP Constructor.
  static SP_base
  Create(SP_input       input,
         SP_mesh        mesh,
         SP_quadrature  quadrature);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set the boundaries for a within-group solve.
  void set(const size_t g);

  /// Update the boundaries for each sweep.
  void update(const size_t g);

  /// brief Update the boundaries for a single angle.
  void update(const size_t g, const size_t o, const size_t a);

  /// Clear the boundary container for a group.
  void clear(const size_t g);

  /// Set the entire group boundary flux for reflecting sides.
  void psi(const size_t g, double *v, const int inout, const int gs,
           bool onlyref = true)
  {
    THROW("IMPLEMENT ME");
  }

  //-------------------------------------------------------------------------//
  // BOUNDARY FLUX ACCESS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Const access to a boundary flux using cardinal indices.
   *
   *  This (and the mutable version) interface is for use
   *  in sweeping, where octants and angles are cycled.
   *
   *  @param    g     Energy group.
   *  @param    o     Octant index.
   *  @param    a     Angle index (within octant).
   *  @param    t     Track index for angle.
   *  @return         Constant reference to boundary flux.
   */
  const double&
  operator()(const size_t g, const size_t o, const size_t a,
             const size_t inout, const size_t t) const;

  /// Mutable access to boundary flux using cardinal indices.
  double&
  operator()(const size_t g, const size_t o, const size_t a,
             const size_t inout, const size_t t);

  /**
   *  @brief Assign the octant, angle, and track fed by an octant,
   *         angle, and track
   *
   *  To make the boundary condition as independent as
   *  possible, BoundaryMOC knows what to do with the
   *  cardinal angle within an octant (i.e. with polar implicit)
   *
   *  @param o1   Incident octant
   *  @param a1   Incident angle within octant
   *  @param t1   Incident track
   *  @param o2   Outgoing octant
   *  @param a2   Outgoing angle within octant
   *  @param t2   Outgoing track
   */
  void feed_into(const size_t o1, const size_t a1, const size_t t1,
                 size_t &o2, size_t &a2, size_t &t2);

  /**
   *  @brief Assign the octant, angle, and track feeding an octant,
   *         angle, and track
   *  @param o1   Outgoing octant
   *  @param a1   Outgoing angle within octant
   *  @param t1   Outgoing track
   *  @param o2   Incident octant
   *  @param a2   Incident angle within octant
   *  @param t2   Incident track
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
  using Base::d_has_vacuum;
  using Base::d_boundary_flux_size;

  /// MOC Quadrature
  SP_quadrature d_quadrature;
  /// Boundary flux [energy, angle, inout, track]
  bf_type d_boundary_flux;
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

#endif // detran_BOUNDARYMOC_HH_

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.hh
//---------------------------------------------------------------------------//
