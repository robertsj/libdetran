//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGPreconditioner.hh
 *  @brief MGPreconditioner
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGPRECONDITIONER_HH_
#define detran_MGPRECONDITIONER_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "transport/FissionSource.hh"
#include "transport/ScatterSource.hh"
#include "utilities/InputDB.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/PCShell.hh"

namespace detran
{

/**
 *  @class MGPreconditioner
 *  @brief Base preconditioner class for multi-group equations
 */

class MGPreconditioner: public callow::PCShell
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef callow::PCShell                           Base;
  typedef detran_utilities::SP<MGPreconditioner>    SP_preconditioner;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef ScatterSource::SP_scattersource           SP_scattersource;
  typedef FissionSource::SP_fissionsource           SP_fissionsource;
  typedef State::SP_state                           SP_state;
  typedef callow::LinearSolver::SP_solver           SP_solver;
  typedef callow::MatrixBase::SP_matrix             SP_operator;
  typedef callow::Vector                            Vector;
  typedef Vector::SP_vector                         SP_vector;
  typedef detran_utilities::size_t                  size_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input            Input database
   *  @param material         Material database
   *  @param mesh             Cartesian mesh
   *  @param ssource          Scattering source
   *  @param fsource          Fission source
   *  @param cutoff           First group included in solve
   *  @param include_fission  Flag to reat fission like scatter
   *  @param adjoint          Flag for adjoint calculations
   *  @param name             Short name for the preconditioner
   */
  MGPreconditioner(SP_input         input,
                   SP_material      material,
                   SP_mesh          mesh,
                   SP_scattersource ssource,
                   SP_fissionsource fsource,
                   size_t           cutoff,
                   bool             include_fission,
                   bool             adjoint,
                   std::string      name = "MG-PC");

  /// virtual destructor
  virtual ~MGPreconditioner(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP PRECONDITIONERS MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  /// Apply the prconditioner to a vector
  virtual void apply(Vector &b, Vector &x) = 0;

  /**
   *  @brief Build the preconditioner
   *
   *  For problems in which materials change (e.g. response generation or
   *  time-dependent models), it may be useful to rebuild the preconditioner
   *  so that it better preconditions the high order operator.
   *
   *  @param  keff    Eigenvalue for scaling fission terms
   */
  virtual void build(const double keff = 1.0, SP_state state = SP_state(0)) = 0;

protected:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Input database
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Cartesian mesh
  SP_mesh d_mesh;
  /// Scatter source
  SP_scattersource d_scattersource;
  /// Fission source
  SP_fissionsource d_fissionsource;
  /// Cutoff
  size_t d_group_cutoff;
  /// Fission flag
  bool d_include_fission;
  /// Adjoint flag
  bool d_adjoint;
  /// Number of groups
  size_t d_number_groups;
  /// Number of active groups
  size_t d_number_active_groups;
  /// Linear solver for inverting the operator
  SP_solver d_solver;
  /// Preconditioning operator
  SP_operator d_operator;
  /// Flag to allow a single build
  bool d_single_build;

};

} // end namespace detran

#endif // detran_MGPRECONDITIONER_HH_

//----------------------------------------------------------------------------//
//              end of file MGPreconditioner.hh
//----------------------------------------------------------------------------//
