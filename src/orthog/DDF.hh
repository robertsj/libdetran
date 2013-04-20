//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   DDF.hh
 *  @brief  DDF
 *  @author Jeremy Roberts
 *  @date   Jan 10, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_orthog_DDF_HH_
#define detran_orthog_DDF_HH_

#include "OrthogonalBasis.hh"

namespace detran_orthog
{

/**
 *  @class DDF
 *  @brief Discrete delta function basis
 *
 *  This initial implementation just piggy-backs on the
 *  default methods.  This could be improved with specialized
 *  fold and transform methods and their inverses.
 */
class DDF: public OrthogonalBasis
{

public:

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *   @brief Constructor.
   *   @param   order   Order of the basis
   *   @param   size    Size of the basis vectors
   */
  DDF(const size_t order, const size_t size)
    : OrthogonalBasis(order, size)
  {
    // Preconditions
    Insist(order + 1 == size,
      "It makes no sense to use DDF for incomplete expansions!");

    // Allocate the basis matrix and fill diagonal.
    d_basis = new callow::MatrixDense(d_size, d_size, 0.0);
    for (size_t i = 0; i < d_order + 1; ++i)
      (*d_basis)(i, i) = 1.0;
  }

  /// Virtual destructor
  virtual ~DDF(){}

};

} // end namespace detran_orthog

#endif // detran_orthog_DDF_HH_

//---------------------------------------------------------------------------//
//              end of file DDF.hh
//---------------------------------------------------------------------------//
