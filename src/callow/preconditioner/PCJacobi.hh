//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PCJacobi.hh
 *  @brief  PCJacobi
 *  @author Jeremy Roberts
 *  @date   Sep 18, 2012
 */
//---------------------------------------------------------------------------//

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

class PCJacobi: public Preconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Preconditioner                   Base;
  typedef typename Base::SP_preconditioner SP_preconditioner;
  typedef typename Matrix::SP_matrix       SP_matrix;
  typedef typename Vector::SP_vector       SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Construct a Jacobi preconditioner for the explicit matrix A
  PCJacobi(SP_matrix A);

  /// Virtual destructor
  virtual ~PCJacobi(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = y
  void apply(Vector &b, Vector &x);

protected:

  /// Inverse diagonal of the user supplied matrix
  SP_vector d_P;

};

} // end namespace callow

#include "PCJacobi.i.hh"

#endif // callow_PCJACOBI_HH_

//---------------------------------------------------------------------------//
//              end of file PCJacobi.hh
//---------------------------------------------------------------------------//
