//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper.hh
 *  @author Jeremy Roberts
 *  @date   Mar 24, 2012
 *  @brief  Sweeper class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_SWEEPER_HH_
#define detran_SWEEPER_HH_

#include "detran_config.hh"
#include "transport/transport_export.hh"
#include "transport/CurrentTally.hh"
#include "transport/Equation.hh"
#include "transport/State.hh"
#include "transport/SweepSource.hh"
#include "angle/Quadrature.hh"
#include "boundary/BoundaryBase.hh"
#include "geometry/Mesh.hh"
#include "material/Material.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"
#include "utilities/SP.hh"
#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Sweeper
 *  @brief Sweeper for discrete ordinates problems.
 *
 *  The within-group transport equation is
 *  @f[
 *      \mathbf{L}\psi = Q \, ,
 *  @f]
 *  where \f$ \mathbf{L} \f$ is the streaming and collision operator and
 *  \f$ Q \f$ is a discrete representation of all source contributions.
 *
 *  To invert the operator \f$ \mathbf{L} \f$, we "sweep" over the mesh for all
 *  angles,
 *  which gives us updated angular fluxes in each cell.  Actually, the
 *  flux *moments* are updated, while the discrete angular flux is
 *  optionally stored.
 *
 *  Relevant input database entries:
 *    - store_angular_flux [int]
 *    - equation [string]
 *
 */
//---------------------------------------------------------------------------//
template <class D>
class Sweeper
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Sweeper>             SP_sweeper;
  typedef State::SP_state                           SP_state;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_angle::Quadrature::SP_quadrature   SP_quadrature;
  typedef BoundaryBase<D>                           Boundary_T;
  typedef typename Boundary_T::SP_boundary          SP_boundary;
  typedef typename SweepSource<D>::SP_sweepsource   SP_sweepsource;
  typedef State::moments_type                       moments_type;
  typedef State::angular_flux_type                  angular_flux_type;
  typedef CurrentTally<D>                           Tally_T;
  typedef typename Tally_T::SP_tally                SP_tally;
  typedef detran_utilities::vec_int                 vec_int;
  typedef detran_utilities::vec2_int                vec2_int;
  typedef detran_utilities::vec3_int                vec3_int;
  typedef detran_utilities::vec_size_t              vec_size_t;
  typedef detran_utilities::vec2_size_t             vec2_size_t;
  typedef detran_utilities::vec3_size_t             vec3_size_t;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  \brief Constructor.
   *
   *  @param    input       User input database.
   *  @param    mesh        Cartesian mesh.
   *  @param    material    Material database.
   *  @param    quadrature  Angular quadrature.
   *  @param    state       State vectors.
   */
  Sweeper(SP_input input,
          SP_mesh mesh,
          SP_material material,
          SP_quadrature quadrature,
          SP_state state,
          SP_boundary boundary,
          SP_sweepsource sweepsource);

  /// Virtual destructor
  virtual ~Sweeper(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SWEEPERS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /**
   *  @brief Sweep over all angles and space.
   *
   *  Note, if the angular flux is to be updated, it is done directly to
   *  via State.  Having sweep take the flux as an explicit argument
   *  allows various input types (e.g. Krylov vectors) without having to go
   *  through State.  The angular flux will never be a direct unknown.
   *
   *  @note __attribute__((always_inline))
   *  
   *  @param phi    reference to moments vector to be updated
   */
  virtual void sweep(moments_type &phi) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Setup the equations for the group.  Default sets only the group index.
  virtual void setup_group(const size_t g);

  /// Allows the psi update to occur whenever needed
  void set_update_psi(const bool v);

  /// Switch on-the-fly boundary updates on or off
  void set_update_boundary(const bool v);

  bool update_boundary() const;

  size_t number_sweeps() const;

  /// Set adjoint
  void set_adjoint(const bool adjoint);

  /// Is adjoint?
  bool is_adjoint() const;

  /// Set a boundary flux tally.
  void set_tally(SP_tally tally);

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Input database
  SP_input d_input;
  /// Mesh
  SP_mesh d_mesh;
  /// Material
  SP_material d_material;
  /// Angular quadrature
  SP_quadrature d_quadrature;
  /// State vectors
  SP_state d_state;
  /// Sweep source
  SP_sweepsource d_sweepsource;
  /// Current group
  size_t d_g;
  /// Update the angular flux?
  bool d_update_psi;
  /// Match incident/outgoing side with octant
  vec3_size_t d_face_index;
  /// Adjoint problem?
  bool d_adjoint;
  /// Count the sweeps.
  size_t d_number_sweeps;
  /// Update the boundary on the fly?  Can't be used for Krylov.
  bool d_update_boundary;
  /// Current tally
  SP_tally d_tally;
  /// Spatial index ranges
  vec3_int d_space_ranges;
  /// Ordered octant indices
  vec_int d_ordered_octants;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /// Allocate template-specific items.
  void setup();

  /// Setup spatial sweep indices.  Avoids functions and handles the adjoint.
  void setup_spatial_indices();

  /// Setup octant sweep indices.
  void setup_octant_indices(SP_boundary);

};

} // end namespace detran

#endif /* detran_SWEEPER_HH_ */
