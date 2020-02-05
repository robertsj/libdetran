//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  WGPreconditioner.hh
 *  @brief WGPreconditioner class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_WGPRECONDITIONER_HH_
#define detran_WGPRECONDITIONER_HH_

#include "material/Material.hh"
#include "geometry/Mesh.hh"
#include "utilities/InputDB.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/PCShell.hh"

namespace detran
{

/**
 *  @class WGPreconditioner
 *  @brief Preconditioner for within-group equation
 *  @note The diffusion operator only operates on the scalar
 *        flux, and so addition of higher order moments will
 *        require restriction and projection operations.
 */

class WGPreconditioner: public callow::PCShell
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef callow::PCShell                           Base;
  typedef detran_utilities::SP<WGPreconditioner>    SP_preconditioner;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_material::Material::SP_material    SP_material;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef callow::LinearSolver::SP_solver           SP_solver;
  typedef std::vector<SP_solver>                    vec_solver;
  typedef callow::MatrixBase::SP_matrix             SP_operator;
  typedef std::vector<SP_operator>                  vec_operator;
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
   */
  WGPreconditioner(SP_input input,
                   SP_material material,
                   SP_mesh mesh,
                   std::string name = "WG-PC");

  /// virtual destructor
  virtual ~WGPreconditioner(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Set the group for this solve.
  void set_group(const size_t group)
  {
    // Preconditions
    Require(group < d_material->number_groups());
    d_group = group;
  }

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
  /// Number of groups
  size_t d_number_groups;
  /// Vector of linear solvers for applying the inverse operator
  vec_solver d_solver;
  /// Vector of one group operators for each group
  vec_operator d_operator;
  /// Group of within group solve
  size_t d_group;

};

} // end namespace detran

#endif // detran_WGPRECONDITIONER_HH_

//---------------------------------------------------------------------------//
//              end of file WGPreconditioner.hh
//---------------------------------------------------------------------------//
