//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  MomentIndexer.hh
 *  @brief MomentIndexer class definition
 *  @note  Copyright (C) Jeremy Roberts 2012-2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_angle_MOMENTINDEXER_HH_
#define detran_angle_MOMENTINDEXER_HH_

#include "angle/angle_export.hh"
#include "utilities/Definitions.hh"
#include "utilities/DBC.hh"
#include "utilities/SP.hh"

namespace detran_angle
{

/**
 *  @class MomentIndexer
 *  @brief Indexex spherical harmonics moments
 *
 *  The one group (or within-group) scattering source in 3-d is defined
 *  @f[
 *      Q(\mathbf{r},\mathbf{\Omega}) =
 *      \sum^L_{l=0} \frac{2l+1}{4\pi}
 *      \sum^l_{m=-l} \Sigma_s{l}(\mathbf{r}) \phi^m_l(\mathbf{r}) \: .
 *  @f]
 *  Such sums over \f$l\f$ and \f$m\f$ are commonplace, in generating
 *  sources and also in simply accessing moments sequentially.  This class
 *  builds a vector of \f$(l,m)\f$ pairs that can be used in a single loop
 *  over all flux moments \f$\phi^m_l\f$.
 *
 *  For 3-d, the moments are ordered as
 *  @f[
 *      [(0,0)] \, , \, \, \, [(1,-1),\, (1,0)\, (1,1)] \, ,
 *      \, \, \, [(2,-2),\,(2,-1)\, \ldots \: .
 *  @f]
 *  For a Legendre order of \f$L\f$, there are \f$(L+1)^2\f$ moments.
 *
 *  For 2-d, several moments can be eliminated by symmetry.  Only those
 *  moments for which \f$l+m\f$ is even are retained (see Hebert).  Thus,
 *  the moments are ordered as
 *  @f[
 *      [(0,0)] \, , \, \, \, [(1,-1),\, (1,1)] \, ,
 *      \, \, \, [(2,-2),\,(2,0)\, \ldots \: .
 *  @f]
 *  For a Legendre order of \f$L\f$, there are \f$(L+1)(L+2)/2\f$ moments.
 *
 *  For 1-d, \f$m=0\f$, so that the moments are ordered by \f$l\f$ alone,
 *  from 0 to the given order.  Hence, there are \f$L+1\f$ moments.
 */

class ANGLE_EXPORT MomentIndexer
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<MomentIndexer>              SP_momentindexer;
  typedef detran_utilities::vec_int                        vec_int;
  typedef detran_utilities::vec2_int                       vec2_int;
  typedef detran_utilities::size_t                         size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param dimension          Problem dimension
   *  @param legendre_order     Legendre order of the \em flux.
   */
  MomentIndexer(const size_t dimension, const size_t legendre_order);

  /// SP constructor
  static SP_momentindexer Create(const size_t dimension,
                                 const size_t legendre_order);

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Return vector of values of m for a given l.
   *  @param  l  Legendre degree.
   *  @return    m values
   */
  const vec_int& m_index(const size_t l) const;

  /// Return legendre order.
  size_t legendre_order() const
  {
    return d_legendre_order;
  }

  /// Return number of moments.
  size_t number_moments() const
  {
    return d_number_moments;
  }

  /**
   *  @brief Return l value.
   *  @param    i   Cardinal moment index.
   *  @return       l value.
   */
  int l(const size_t i) const;

  /**
   *  @brief Return m value.
   *  @param    i   Cardinal moment index.
   *  @return       m value.
   */
  int m(const size_t i) const;

  /**
   *  @brief Return l value.
   *  @param    l   Cardinal moment index.
   *  @param    m   Cardinal moment index.
   *  @return       Cardinal moment index.
   */
  size_t index(const size_t l, const int m) const;

  /// Print the indices
  void display() const;

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Legendre order of the \em flux.
  size_t d_legendre_order;
  /// Number of moments.
  size_t d_number_moments;
  /// Vector of proper \f$m\f$ values for a given \f$l\f$.
  vec2_int d_m_index;
  /// Vector of l values.
  vec_int d_l;
  /// Vector of m values.
  vec_int d_m;

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  void construct_1D();
  void construct_2D();
  void construct_3D();

};

ANGLE_TEMPLATE_EXPORT(detran_utilities::SP<MomentIndexer>)

} // end namespace detran_angle

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "angle/MomentIndexer.i.hh"


#endif /* MOMENTINDEXER_HH_ */
