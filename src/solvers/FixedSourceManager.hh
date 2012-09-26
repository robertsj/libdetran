//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FixedSourceManager.hh
 *  @author robertsj
 *  @date   Sep 25, 2012
 *  @brief  FixedSourceManager class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_FIXEDSOURCEMANAGER_HH_
#define detran_FIXEDSOURCEMANAGER_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "callow/solver/LinearSolver.hh"
#include "transport/State.hh"
#include "utilities/Definitions.hh"
#include "boundary/BoundaryBase.hh"
#include <vector>

namespace detran
{

/**
 *  @class FixedSourceManager
 *  @brief Manage solution of a multigroup fixed source problem
 */
template <class D>
class FixedSourceManager
{

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::InputDB::SP_input         SP_input;
  typedef State::SP_state                             SP_state;
  typedef detran_geometry::Mesh::SP_mesh              SP_mesh;
  typedef detran_material::Material::SP_material      SP_material;
  //
  typedef BoundaryBase<D>                             Boundary_T;
  typedef typename Boundary_T::SP_boundary            SP_boundary;
  //
  typedef detran_external_source::
          ExternalSource::SP_externalsource           SP_source;
  //
  typedef State::moments_type                         moments_type;
  //
  typedef callow::LinearSolver<double>::SP_solver     SP_solver;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      parameter database
   *  @param material   material database
   *  @param mesh       mesh definition
   *  @param state      state vector
   */
  FixedSourceManager(SP_input    input,
                     SP_material material,
                     SP_mesh     mesh,
                     SP_state    state);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Add a new external source
  void set_source(SP_source q);

  /// Solve the system
  void solve();

  /// @name Getters
  /// @{
  SP_input input() const { return d_input; }
  SP_material material() const { return d_material; }
  SP_mesh mesh() const { return d_mesh; }
  SP_state state() const { return d_state; }
  SP_boundary boundary() const { return d_boundary; }
  /// @}

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// User input
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Mesh
  SP_mesh d_mesh;
  /// State vector
  SP_state d_state;
  /// Boundary container
  SP_boundary d_boundary;
  /// External volume sources
  std::vector<SP_source> d_sources;
  /// External boundary sources
  //std::vector<SP_source> d_sources;

  /// Linear solver
  SP_solver d_solver;

  /// Adjoint mode flag
  bool d_adjoint;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//


};

} // end namespace detran

#endif /* detran_FIXEDSOURCEMANAGER_HH_ */
