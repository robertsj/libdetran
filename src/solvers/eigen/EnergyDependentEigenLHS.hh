//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EnergyDependentEigenLHS.hh
 *  @brief EnergyDependentEigenLHS class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_ENERGYDEPENDENTEIGENLHS_HH_
#define detran_ENERGYDEPENDENTEIGENLHS_HH_

#include "solvers/FixedSourceManager.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class EnergyDependentEigenLHS
 *  @brief Right hand side of energy-dependent eigenvalue problem
 *
 *  This operator represents the action of
 *  @f[
 *      \mathbf{TMS}\boldsymbol{\chi}\mathbf{F}^T \, .
 *  @f]
 *
 *  The left hand operator is represented by MGTransportOperator.
 */

template <class D>
class EnergyDependentEigenLHS: public callow::MatrixShell
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MatrixShell                                   Base;
  typedef EnergyDependentEigenLHS<D>                    Operator_T;
  typedef typename detran_utilities::SP<Operator_T>     SP_operator;
  typedef FixedSourceManager<D>                         Fixed_T;
  typedef typename Fixed_T::SP_manager                  SP_mg_solver;
  typedef typename Fixed_T::SP_state                    SP_state;
  typedef typename Fixed_T::SP_fissionsource            SP_fissionsource;
  typedef callow::Vector                                Vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param mg_solver  Multigroup solver
   */
  EnergyDependentEigenLHS(SP_mg_solver mg_solver);

  // Destructor
  virtual ~EnergyDependentEigenLHS(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  // Default shell display gives just the sizes
  virtual void display() const;

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  virtual void multiply(const Vector &x,  Vector &y);

  // the client must implement the action y <-- A' * x
  virtual void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("Transpose transport operator not implemented");
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Multigroup solver
  SP_mg_solver d_mg_solver;
  /// Fission source
  SP_fissionsource d_fissionsource;

};

} // end namespace detran




#endif /* detran_ENERGYDEPENDENTEIGENLHS_HH_ */
