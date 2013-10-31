//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  MGTransportOperator.hh
 *  @brief MGTransportOperator class definition
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef detran_MGTRANSPORTOPERATOR_HH_
#define detran_MGTRANSPORTOPERATOR_HH_

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
 *  @class MGTransportOperator
 *  @brief Multigroup transport operator
 *
 *  The multigroup transport operator is defined as
 *  ... finish.
 */

template <class D>
class MGTransportOperator: public callow::MatrixShell
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef MatrixShell                                 Base;
  typedef MGTransportOperator<D>                      Operator_T;
  typedef typename detran_utilities::SP<Operator_T>   SP_operator;
  typedef State::SP_state                             SP_state;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_utilities::size_t                    size_t;
  typedef detran_utilities::vec_size_t                groups_t;
  typedef groups_t::iterator                          groups_iter;
  typedef typename Sweeper<D>::SP_sweeper             SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;
  typedef State::moments_type                         moments_type;
  typedef callow::Vector                              Vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state      state vector
   *  @param boundary   boundary container
   *  @param sweeper    transport sweeper
   *  @param source     sweep source
   *  @param cutoff     lowest group included in operator
   *  @param adjoint    adjoint flag
   */
  MGTransportOperator(SP_state        state,
                      SP_boundary     boundary,
                      SP_sweeper      sweeper,
                      SP_sweepsource  source,
                      size_t          cutoff,
                      bool            adjoint);

  // Destructor
  virtual ~MGTransportOperator(){}

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  // default shell display gives just the sizes
  virtual void display() const;

  size_t moments_size() const { return d_moments_size; }
  size_t boundary_size() const { return d_boundary_size; }

  SP_state state() { return d_state; }
  SP_boundary boundary() { return d_boundary; }
  SP_sweeper sweeper() { return d_sweeper; }
  SP_sweepsource sweepsource() { return d_sweepsource; }


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

#endif // detran_MGTRANSPORTOPERATOR_HH_

//----------------------------------------------------------------------------//
//              end of file MGTransportOperator.hh
//----------------------------------------------------------------------------//
