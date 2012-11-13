//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGPreconditioner.hh
 *  @brief  MGPreconditioner
 *  @author Jeremy Roberts
 *  @date   Nov 12, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGPRECONDITIONER_HH_
#define detran_MGPRECONDITIONER_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/PCShell.hh"

namespace detran
{

/**
 *  @class MGPreconditioner
 *  @brief Base preconditioner class for multi-group equations
 *  @note The diffusion operator only operates on the scalar
 *        flux, and so addition of higher order moments will
 *        require restriction and projection operations.
 */

class MGPreconditioner: public callow::PCShell
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef callow::PCShell                           Base;
  typedef detran_utilities::SP<MGPreconditioner>    SP_preconditioner;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef callow::LinearSolver::SP_solver           SP_solver;
  typedef callow::MatrixBase::SP_matrix             SP_operator;
  typedef callow::Vector                            Vector;
  typedef Vector::SP_vector                         SP_vector;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input      Input database
   *  @param material   Material database
   *  @param mesh       Cartesian mesh
   *  @param cutoff     Lowest group to include in the operator
   */
  MGPreconditioner(SP_input input,
                   SP_material material,
                   SP_mesh mesh,
                   size_t cutoff,
                   std::string name = "MG-PC");

  /// virtual destructor
  virtual ~MGPreconditioner(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  virtual void apply(Vector &b, Vector &x) = 0;

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Input database
  SP_input d_input;
  /// Material database
  SP_material d_material;
  /// Cartesian mesh
  SP_mesh d_mesh;
  /// Cutoff
  size_t d_group_cutoff;
  /// Number of groups
  size_t d_number_groups;
  /// Number of active groups
  size_t d_number_active_groups;
  /// Vector of linear solvers for applying the inverse operator
  SP_solver d_solver;
  /// Vector of diffusion loss operators for each group
  SP_operator d_operator;

};

} // end namespace detran

#endif // detran_MGPRECONDITIONER_HH_

//---------------------------------------------------------------------------//
//              end of file MGPreconditioner.hh
//---------------------------------------------------------------------------//
