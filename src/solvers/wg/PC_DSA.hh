//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PC_DSA.hh
 *  @brief PC_DSA class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_PC_DSA_HH_
#define detran_PC_DSA_HH_

#include "WGPreconditioner.hh"
#include "WGDiffusionLossOperator.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"

namespace detran
{

/**
 *  @class PC_DSA
 *  @brief Diffusion synthetic acceleration
 *
 *  The DSA preconditioning process \$ \mathbf{P}^{-1} \$
 *  is defined to be
 *  @f[
 *      (\mathbf{I} - \mathbf{C}^{-1} \mathbf{S}) \, ,
 *  @f]
 *  where \f$ \mathbf{C} \f$ is the one group diffusion operator.
 *
 *  @todo Include fission if treated like scatter
 */

class PC_DSA: public WGPreconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef WGPreconditioner                  Base;
  typedef detran_utilities::SP<PC_DSA>      SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef WGDiffusionLossOperator           Operator_T;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

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
   */
  PC_DSA(SP_input input,
         SP_material material,
         SP_mesh mesh,
         SP_scattersource source);

  /// virtual destructor
  virtual ~PC_DSA(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// solve Px = b
  void apply(Vector &b, Vector &x);

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

//----------------------------------------------------------------------------//
// INLINE FUNCTIONS
//----------------------------------------------------------------------------//

#include "PC_DSA.i.hh"

#endif // detran_PC_DSA_HH_

//----------------------------------------------------------------------------//
//              end of file PC_DSA.hh
//----------------------------------------------------------------------------//
