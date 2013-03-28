//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   CMMGDSA.hh
 *  @brief  CMMGDSA 
 *  @author Jeremy Roberts
 *  @date   Mar 26, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_CMMGDSA_HH_
#define detran_CMMGDSA_HH_

#include "MGPreconditioner.hh"
#include "DiffusionLossOperator.hh"
#include "transport/CoarseMesh.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class CMMGDSA
 *  @brief Coarse mesh, multigroup diffusion synthetic acceleration
 *
 *  The multigroup DSA preconditioning process \$ \mathbf{P}^{-1} \$
 *  is defined to be
 *  @f[
 *      (\mathbf{I} - \mathbf{R}^T \mathbf{C}^{-1} \mathbf{R} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the multigroup diffusion operator on
 *  a coarse spatial mesh and
 *  \f$ \mathbf{R} \f$ and its transpose represent a spatial restriction
 *  and projection, respectively.
 *  This operator treats group-to-group scattering and, if requested,
 *  fission implicitly.
 *
 *
 *  @note This inherits from the shell matrix (for now) so that the
 *        action can be used to construct an explicit operator for
 *        detailed numerical studies
 */

class CMMGDSA: public callow::MatrixShell, public MGPreconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<CMMGDSA>     SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef DiffusionLossOperator             Operator_T;
  typedef CoarseMesh::SP_coarsemesh         SP_coarsemesh;
  typedef CoarseMesh::SP_mesh               SP_mesh;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  Assuming the within-group transport problem is set up,
   *  a KSP object exists from which the PC is extracted.  This
   *  PC is passed here to be constructed and for its application
   *  operator to be assigned.
   *
   *  @param input      Input database
   *  @param material   Material database
   *  @param mesh       Cartesian mesh
   *  @param source     Scattering source
   *  @param cutoff     First group included in solve
   */
  CMMGDSA(SP_input          input,
          SP_material       material,
          SP_mesh           mesh,
          SP_scattersource  source,
          size_t            cutoff,
          bool              include_fission);

  /// virtual destructor
  virtual ~CMMGDSA(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// solve Px = b
  void apply(Vector &b, Vector &x);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SHELL MATRICES MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  void multiply(const Vector &x,  Vector &y)
  {
    Vector b(x.size(), 0.0);
    b.copy(x);
    apply(b, y);
  }

  // the client must implement the action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("NOT IMPLEMENTED");
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Coarse mesh


  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

#endif // detran_CMMGDSA_HH_

//---------------------------------------------------------------------------//
//              end of file CMMGDSA.hh
//---------------------------------------------------------------------------//
