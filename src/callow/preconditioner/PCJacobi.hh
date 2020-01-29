//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCJacobi.hh
 *  @brief PCJacobi class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_PCJACOBI_HH_
#define callow_PCJACOBI_HH_

#include "Preconditioner.hh"
#include "callow/matrix/Matrix.hh"

namespace callow
{

/**
 *  @class PCJacobi
 *  @brief Applies a Jacobi preconditioner
 *
 *  The Jacobi preconditioner is defined by the process
 *  @f[
 *      \mathbf{P}^{-1} = \mathbf{D}^{-1}
 *                      = \mathrm{diag}([a^{-1}_{11}, a^{-1}_{22}, \cdots]^T) \, ,
 *  @f]
 *  where \f$ a_{ii} \f$ is the \e ith diagonal element of \f$ \mathbf{A} \f$.
 *
 *  If zero elements are found on the diagonal, the inverse is set to the
 *  matrix size.
 *
 */

class CALLOW_EXPORT PCJacobi: public Preconditioner
{

public:

  //--------------------------------------------------------------------------//
  // TYPEDEFS
  //--------------------------------------------------------------------------//

  typedef Preconditioner                    Base;
  typedef Base::SP_preconditioner           SP_preconditioner;
  typedef MatrixBase::SP_matrix             SP_matrix;
  typedef Matrix::SP_matrix                 SP_matrixfull;
  typedef Vector::SP_vector                 SP_vector;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  /// Construct a Jacobi preconditioner for the explicit matrix A
  PCJacobi(SP_matrix A);

  /// SP constructor
  static SP_preconditioner Create(SP_matrix A);

  /// Virtual destructor
  virtual ~PCJacobi(){};

  //--------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //--------------------------------------------------------------------------//

  /// Solve Px = y
  void apply(Vector &b, Vector &x);

protected:

  /// Inverse diagonal of the user supplied matrix
  SP_vector d_P;

};

} // end namespace callow

#include "PCJacobi.i.hh"

#endif // callow_PCJACOBI_HH_

//----------------------------------------------------------------------------//
//              end of file PCJacobi.hh
//----------------------------------------------------------------------------//
