//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundarySN.hh
 *  @author Jeremy Roberts
 *  @date   Mar 24, 2012
 *  @brief  BoundarySN class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYSN_HH_
#define detran_BOUNDARYSN_HH_

#include "BoundaryBase.hh"
#include "BoundaryTraits.hh"
#include "BoundaryCondition.hh"
#include "angle/Quadrature.hh"

namespace detran
{

template <class D> class BoundaryCondition;

//---------------------------------------------------------------------------//
/**
 *  @class BoundarySN
 *  @brief Boundary flux container for SN problems.
 *  @todo Switch accessor to match diffusion (or vv)
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
  typedef typename BoundaryTraits<D>::value_type    bf_type;
  typedef typename std::vector<bf_type>             vec_boundary_flux;
  typedef typename std::vector<vec_boundary_flux>   vec2_boundary_flux;
  typedef typename std::vector<vec2_boundary_flux>  vec3_boundary_flux;
  typedef typename BoundaryCondition<D>::SP_bc      SP_bc;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_utilities::vec_dbl                 vec_dbl;
  typedef detran_utilities::size_t                  size_t;
  typedef D                                         D_T;

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
   *  @param    quadrature  Angular quadrature.
   */
  BoundarySN(SP_input        input,
             SP_mesh         mesh,
             SP_quadrature   quadrature);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set the boundaries for a within-group solve.
  void set(const size_t g);

  /// Update the boundaries for each sweep.
  void update(const size_t g);

  /// Update the boundaries for a single angle.
  void update(const size_t g, const size_t o, const size_t a);

  /// Clear the boundary container for a group.
  void clear(const size_t g);

  /// Set the entire group boundary flux for reflecting sides.
  void psi(const size_t g, double *v, const int inout, const int gs,
           bool onlyref = true);

  //-------------------------------------------------------------------------//
  // BOUNDARY FLUX ACCESS
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
  operator()(const size_t side,
             const size_t o,
             const size_t a,
             const size_t g) const;

  /// Mutable access to boundary flux using cardinal indices.
  bf_type&
  operator()(const size_t side,
             const size_t o,
             const size_t a,
             const size_t g);

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

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Clear any fixed boundary conditions.
  void clear_bc()
  {
    for (size_t s = 0; s < d_bc.size(); ++s)
    {
      d_bc[s]->clear();
    }
  }

  /// Display boundray information and contents
  void display(bool inout) const;

  /// Set the boundary condition
  void set_bc(const size_t side, SP_bc b)
  {
    Require(side < D::dimension*2);
    Require(b);
    d_bc[side] = b;
  }

  /// Get the boundary condition.
  SP_bc bc(const size_t side)
  {
    Require(side < D::dimension*2);
    return d_bc[side];
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

#endif /* detran_BOUNDARYSN_HH_ */

//---------------------------------------------------------------------------//
//              end of BoundarySN.hh
//---------------------------------------------------------------------------//
