//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  SweepOperator.hh
 *  @brief SweepOperator
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_SWEEPOPERATOR_HH_
#define detran_SWEEPOPERATOR_HH_

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
 *  @class SweepOperator
 *  @brief Operator wrapped around a transport sweep
 */
template <class D>
class SweepOperator: public callow::MatrixShell
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MatrixShell                                 Base;
  typedef SweepOperator<D>                            Operator_T;
  typedef typename detran_utilities::SP<Operator_T>   SP_operator;
  typedef State::SP_state                             SP_state;
  typedef typename BoundaryBase<D>::SP_boundary       SP_boundary;
  typedef detran_utilities::size_t                    size_t;
  typedef typename Sweeper<D>::SP_sweeper             SP_sweeper;
  typedef typename SweepSource<D>::SP_sweepsource     SP_sweepsource;
  typedef State::moments_type                         moments_type;
  typedef callow::Vector                              Vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  // Constructor
  SweepOperator(SP_state        state,
                SP_boundary     boundary,
                SP_sweeper      sweeper,
                SP_sweepsource  source);

  // Destructor
  virtual ~SweepOperator(){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  // default shell display gives just the sizes
  virtual void display() const;

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MATRICES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  // the client must implement the action y <-- A * x
  virtual void multiply(const Vector &x,  Vector &y);

  // the client must implement the action y <-- A' * x
  virtual void multiply_transpose(const Vector &x, Vector &y)
  {
    THROW("Transpose sweep operator not implemented");
  }

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  SP_state d_state;
  SP_boundary d_boundary;
  SP_sweeper d_sweeper;
  SP_sweepsource d_sweepsource;
  size_t d_moments_size;
  size_t d_boundary_size;
  size_t d_g;

};


} // end namespace detran

#endif // detran_SWEEPOPERATOR_HH_

//---------------------------------------------------------------------------//
//              end of file SweepOperator.hh
//---------------------------------------------------------------------------//
