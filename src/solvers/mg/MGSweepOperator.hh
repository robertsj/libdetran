//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGSweepOperator.hh
 *  @brief MGSweepOperator class definition
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//


#ifndef detran_MGSweepOperator_HH_
#define detran_MGSweepOperator_HH_

#include "boundary/BoundaryBase.hh"
#include "callow/matrix/MatrixShell.hh"
#include "transport/State.hh"
#include "transport/Sweeper.hh"
#include "transport/SweepSource.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"

namespace detran
{

/**
 *  @class MGSweepOperator
 *  @brief Performs one multigroup sweep, denoted
 *         @f$ \mathbf{DL}^{-1} \mathbf{M} @f$.
 */

template <class D>
class MGSweepOperator: public callow::MatrixShell
{

public:


  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MatrixShell                                 Base;
  typedef MGSweepOperator                             Operator_T;
  typedef detran_utilities::SP<Operator_T>            SP_operator;
  typedef State::SP_state                             SP_state;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_utilities::size_t                    size_t;
  typedef detran_utilities::vec_size_t                groups_t;
  typedef groups_t::iterator                          groups_iter;
  typedef typename Sweeper<D>::SP_sweeper             SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;
  typedef State::moments_type                         moments_type;
  typedef callow::Vector                              Vector;

  /**
   *  @brief Constructor
   */
  MGSweepOperator(SP_state        state,
                  SP_boundary     boundary,
                  SP_sweeper      sweeper,
                  SP_sweepsource  source,
                  size_t          cutoff,
                  bool            adjoint);

  /// Virtual destructor
  virtual ~MGSweepOperator(){}

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT THESE
  //--------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  void multiply(const Vector &x,  Vector &y);

  // the client must implement the action y <-- A' * x
  void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("Sweep transpose operator not implemented");
  }

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// State vector
  SP_state d_state;
  /// Boundary flux container
  SP_boundary d_boundary;
  /// Transport sweeper
  SP_sweeper d_sweeper;
  /// Sweep source
  SP_sweepsource d_sweepsource;
  /// Total number of groups
  size_t d_number_groups;
  /// Number of groups for which operator is applicable
  size_t d_number_active_groups;
  /// Group index below which the operator is not applicable
  size_t d_krylov_group_cutoff;
  /// Size of a group moment vector
  size_t d_moments_size;
  /// Size of a group boundary vector
  size_t d_boundary_size;
  /// Adjoint flag
  bool d_adjoint;
  /// Lower group bound
  groups_t d_groups;


};

} // end namespace detran

#endif /* detran_MGSweepOperator_HH_ */

//----------------------------------------------------------------------------//
//              end of file MGSweepOperator.hh
//----------------------------------------------------------------------------//
