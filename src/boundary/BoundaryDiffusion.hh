//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryDiffusion.hh
 *  @brief  BoundaryDiffusion class definition
 *  @author Jeremy Roberts
 *  @date   Sep 10, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYDIFFUSION_HH_
#define detran_BOUNDARYDIFFUSION_HH_

#include "BoundaryBase.hh"
#include "BoundaryTraits.hh"

namespace detran
{

/**
 *  @class BoundaryDiffusion
 *  @brief Container for diffusion boundary partial currents
 *
 *  Unlike the boundary containers for transport, the diffusion boundary
 *  is not responsible for boundary conditions.  This is because these
 *  conditions are built into the diffusion operator.
 */

template <class D>
class BoundaryDiffusion: public BoundaryBase<D>
{
public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryBase<D>                               Base;
  typedef typename Base::SP_boundary                    SP_base;
  typedef detran_utilities::SP<BoundaryDiffusion<D> >   SP_boundary;
  typedef typename Base::SP_input                       SP_input;
  typedef typename Base::SP_mesh                        SP_mesh;
  typedef typename BoundaryTraits<D>::value_type        bf_type;
  typedef typename std::vector<bf_type>                 vec_boundary_flux;
  typedef typename std::vector<vec_boundary_flux>       vec2_boundary_flux;
  typedef typename std::vector<vec2_boundary_flux>      vec3_boundary_flux;
  typedef detran_geometry::Mesh                         Mesh;
  typedef detran_utilities::vec_dbl                     vec_dbl;
  typedef detran_utilities::size_t                      size_t;
  typedef D                                             D_T;

  using Base::IN;
  using Base::OUT;

  //-------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *
   *  @param    input       User input database.
   *  @param    mesh        Cartesian mesh.
   */
  BoundaryDiffusion(SP_input        input,
                    SP_mesh         mesh);

  /// SP Constructor
  static SP_base
  Create(SP_input         input,
         SP_mesh          mesh);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  @param  g   Group of current solve
   */
  void set(const size_t g){ /* not needed */ }

  /**
   *  @brief Update the boundaries for each sweep.
   *
   *  This updates all incident boundary fluxes
   *  using the current outgoing boundary fluxes
   *  in a group.  What happens in the update is
   *  a function of the underlying boundary
   *  condition.
   *
   *  @param  g   Group of current solve
   */
  void update(const size_t g){ /* not needed */ }

  /**
   *  @brief Update the boundaries for a single angle.
   *
   *  This is an alternative update that only updates
   *  the incident boundary flux for a particular
   *  angle.  When called within a sweep, this allows
   *  the most recent boundary fluxes to be used,
   *  in effect producing a Gauss-Seidel iteration.
   *
   *  @note This cannot be used for Krylov solvers.
   *
   *  @param  g   Group of current solve
   *  @param  o   Octant being swept
   *  @param  a   Angle within octant being swept
   */
  void update(const size_t g, const size_t o, const size_t a)
  { /* not needed */ }

  /**
   *  @brief Clear the boundary container for a group.
   *
   *  In some cases, a client might require homogeneous
   *  boundaries, perhaps after a fixed boundary has
   *  been used to construct a right hand side for a
   *  Krylov solve.
   *
   *  @param  g   Group of current solve
   */
  void clear(const size_t g);

  /**
   *  @brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is is support of Krylov solvers.
   *
   *  @param g  Group of current sweep
   *  @param v  Pointer to data used in Krylov solve
   */
  void psi(const size_t g, double *v, const int inout, const int gs,
           bool onlyref = true)
  { /* not needed */ }


  //-------------------------------------------------------------------------//
  // BOUNDARY FLUX ACCES
  //-------------------------------------------------------------------------//

  /**
   *  @brief Const access to a boundary flux using cardinal indices.
   *
   *  This (and the mutable version) interface is for use
   *  in sweeping, where octants and angles are cycled.
   *
   *  @param    side  Side index.
   *  @param    o     Octant index.
   *  @param    a     Angle index (within octant).
   *  @param    g     Energy group.
   *  @return         Constant reference to boundary flux.
   */
  const bf_type&
  operator()(const size_t side, const size_t g, const size_t inout) const;

  /// Mutable access to boundary flux using cardinal indices.
  bf_type&
  operator()(const size_t side, const size_t g, const size_t inout);

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

  /// Boundary flux (inout, side, energy).(space^D)
  vec3_boundary_flux d_boundary_flux;

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

#include "BoundaryDiffusion.i.hh"

#endif // detran_BOUNDARYDIFFUSION_HH_

//---------------------------------------------------------------------------//
//              end of file BoundaryDiffusion.hh
//---------------------------------------------------------------------------//
