//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MGCMDSA.hh
 *  @brief MGCMDSA class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGCMDSA_HH_
#define detran_MGCMDSA_HH_

#include "MGCoarseMeshPreconditioner.hh"
#include "DiffusionLossOperator.hh"

namespace detran
{

/**
 *  @class MGCMDSA
 *  @brief Coarse mesh, multigroup diffusion synthetic acceleration
 *
 *  The multigroup DSA preconditioning process \$ \mathbf{P}^{-1} \$
 *  is defined to be
 *  @f[
 *      (\mathbf{I} + \mathbf{P} \mathbf{C}^{-1} \mathbf{R} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the multigroup diffusion operator on
 *  a coarse spatial mesh and
 *  \f$ \mathbf{R} \f$ and  \f$ \mathbf{P} \f$ represent the space-energy
 *  restriction and prolongation operators, respectively.
 *
 *  See @ref MGCoarseMeshPreconditioner for more details on defining the
 *  coarse space and energy meshes.
 */

class MGCMDSA: public MGCoarseMeshPreconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGCoarseMeshPreconditioner        Base;
  typedef DiffusionLossOperator             Operator_T;
  typedef callow::Matrix                    Matrix;
  typedef Matrix::SP_matrix                 SP_matrix;

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
          SP_scattersource  ssource,
          SP_fissionsource  fsource,
          size_t            cutoff,
          bool              include_fission);

  /// virtual destructor
  virtual ~MGCMDSA(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP PRECONDITIONERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// solve Px = b
  void apply(Vector &b, Vector &x);

  /// build the preconditioner
  void build(const double keff = 1.0, SP_state state = SP_state(0));


private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Coarse mesh
  SP_coarsemesh d_coarsemesher;
  /// Scatter source
  SP_scattersource d_scattersource;
  /// Restriction operator
  SP_matrix d_restrict;
  /// Projection operator
  SP_matrix d_project;
  /// Flag to include fission
  bool d_include_fission;

};

} // end namespace detran

#endif // detran_MGCMDSA_HH_

//---------------------------------------------------------------------------//
//              end of file MGCMDSA.hh
//---------------------------------------------------------------------------//
