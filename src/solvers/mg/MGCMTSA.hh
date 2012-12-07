//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   MGCMTSA.hh
 *  @author robertsj
 *  @date   Dec 6, 2012
 *  @brief  MGCMTSA class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGCMTSA_HH_
#define detran_MGCMTSA_HH_

#include "MGPreconditioner.hh"
#include "DiffusionLossOperator.hh"
#include "transport/ScatterSource.hh"
#include "callow/solver/LinearSolver.hh"
#include "callow/preconditioner/Preconditioner.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class MGCMTSA.hh
 *  @brief Multigroup coarse mesh transport synthetic acceleration
 */

class MGCMTSA: public MGPreconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGPreconditioner                  Base;
  typedef detran_utilities::SP<MGCMTSA>     SP_pc;
  typedef ScatterSource::SP_scattersource   SP_scattersource;
  typedef DiffusionLossOperator             Operator_T;

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
  MGCMTSA(SP_input input,
        SP_material material,
        SP_mesh mesh,
        SP_scattersource source,
        size_t cutoff,
        bool include_fission);

  /// virtual destructor
  virtual ~MGCMTSA(){}

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

  /// Scatter source
  SP_scattersource d_scattersource;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//


#endif /* MGCMTSA_HH_ */
