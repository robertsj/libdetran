//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MGCMDSA.hh
 *  @brief MGCMDSA class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGCMDSA_HH_
#define detran_MGCMDSA_HH_

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
 *  @class MGCMDSA
 *  @brief Coarse mesh, multigroup diffusion synthetic acceleration
 *
 *  The multigroup DSA preconditioning process \$ \mathbf{P}^{-1} \$
 *  is defined to be
 *  @f[
 *      (\mathbf{I} - \mathbf{P} \mathbf{C}^{-1} \mathbf{R} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the multigroup diffusion operator on
 *  a coarse spatial mesh and
 *  \f$ \mathbf{R} \f$ and  \f$ \mathbf{P} \f$ represent the space-energy
 *  restriction and projection operators, respectively.
 *
 *  The restriction space is based on a user-specific level defining the
 *  number of fine cells per coarse cell.  The same is true for energy.
 *  The current implementation uses only unity-weighted cross section
 *  collapsing.  It would be possible to use some sort of approximate
 *  spectral shape as well based on infinite medium calculations in one
 *  or more material mixture.
 *
 *  The restriction in space is based either on a user-defined state vector,
 *  possibly from an initial guess, partial solution, or other approximation,
 *  or uniform volume weighting.
 *
 *  The restriction in energy is based either on the same user defined state
 *  vector or a set of material-dependent spectra based on infinite medium
 *  calculations for each material.
 */

class MGCMDSA: public MGPreconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<MGCMDSA>     SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef DiffusionLossOperator             Operator_T;
  typedef CoarseMesh::SP_coarsemesh         SP_coarsemesh;
  typedef CoarseMesh::SP_mesh               SP_mesh;
  typedef State::SP_state                   SP_state;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *
   *  @param input            Input database
   *  @param material         Material database
   *  @param mesh             Cartesian mesh
   *  @param source           Scattering source
   *  @param cutoff           First group included in solve
   *  @param include_fission  Treat fission like scatter
   */
  MGCMDSA(SP_input          input,
          SP_material       material,
          SP_mesh           mesh,
          SP_scattersource  source,
          size_t            cutoff,
          bool              include_fission);

  /// virtual destructor
  virtual ~MGCMDSA(){}

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
  SP_coarsemesh d_coarsemesher;

  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

#endif // detran_MGCMDSA_HH_

//---------------------------------------------------------------------------//
//              end of file MGCMDSA.hh
//---------------------------------------------------------------------------//
