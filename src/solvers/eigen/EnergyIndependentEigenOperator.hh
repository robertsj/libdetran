//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  EnergyIndependentEigenOperator.hh
 *  @brief EnergyIndependentEigenOperator class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_ENERGYINDEPENDENTEIGENOPERATOR_HH_
#define detran_ENERGYINDEPENDENTEIGENOPERATOR_HH_

#include "solvers/FixedSourceManager.hh"
#include "callow/matrix/MatrixShell.hh"

namespace detran
{

/**
 *  @class EnergyIndependentEigenOperator
 *  @brief Energy-independent operator for eigenvalue problems
 *
 *  This operator represents the action of
 *  @f[
 *      \mathbf{F}^T (\mathbf{I} - \mathbf{TMS})^{-1}
 *      \mathbf{TM} \boldsymbol{\chi} \, .
 *  @f]
 *
 *  This is the preferred operator for transport eigenvalue
 *  solves since the system is reduced in size.
 *
 */
template <class D>
class EnergyIndependentEigenOperator: public callow::MatrixShell
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MatrixShell                                   Base;
  typedef EnergyIndependentEigenOperator<D>             Operator_T;
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
  EnergyIndependentEigenOperator(SP_mg_solver mg_solver);

  // Destructor
  virtual ~EnergyIndependentEigenOperator(){}

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

#endif /* detran_ENERGYINDEPENDENTEIGENOPERATOR_HH_ */

//----------------------------------------------------------------------------//
//              end of file EnergyIndependentEigenOperator.hh
//----------------------------------------------------------------------------//
