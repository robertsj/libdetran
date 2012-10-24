//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   PCIdentity.hh
 *  @author robertsj
 *  @date   Oct 24, 2012
 *  @brief  PCIdentity class definition.
 */
//---------------------------------------------------------------------------//

#ifndef callow_PCIDENTITY_HH_
#define callow_PCIDENTITY_HH_

#include "Preconditioner.hh"

namespace callow
{

/**
 *  @class PCIdentity
 *  @brief Implements an indentity preconditioner (i.e. no preconditioning)
 *
 *  This is mostly for testing purposes.
 */
class PCIdentity: public Preconditioner
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Preconditioner                    Base;
  typedef Base::SP_preconditioner           SP_preconditioner;
  typedef MatrixBase::SP_matrix             SP_matrix;
  typedef Matrix::SP_matrix                 SP_matrixfull;
  typedef Vector::SP_vector                 SP_vector;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  PCIdentity(double factor = 1.0)
    : d_factor(factor)
  {/* ... */}

  /// SP constructor
  static SP_preconditioner Create(double factor = 1.0)
  {
    SP_preconditioner p(new PCIdentity(factor));
    return p;
  }

  /// Virtual destructor
  virtual ~PCIdentity(){};

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL PRECONDITIONERS MUST IMPLEMENT THIS
  //-------------------------------------------------------------------------//

  /// Solve Px = b
  void apply(Vector &b, Vector &x)
  {
    x.copy(b);
    x.scale(d_factor);
  }

private:

  /// Scaling factor
  double d_factor;

};

} // end namespace callow

#endif /* callow_PCIDENTITY_HH_ */
