//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussLegendre.hh
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 * \brief  GaussLegendre class definition.
 */
//---------------------------------------------------------------------------//

#ifndef GAUSSLEGENDRE_HH_
#define GAUSSLEGENDRE_HH_

// Angle headers
#include "Quadrature.hh"

namespace detran_angle
{

class GaussLegendre : public Quadrature
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<GaussLegendre> SP_quadrature;
  typedef detran_utilities::SP<Quadrature>    SP_base;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param    order       Quadrature order.
   */
  GaussLegendre(size_t order);

  /// SP constructor
  static detran_utilities::SP<Quadrature> Create(size_t order)
  {
    SP_quadrature p(new GaussLegendre(order));
    return p;
  }

private:

  //-------------------------------------------------------------------------//
  // IMPLEMENTATION
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Compute GL parameters for arbitrary order.
   *  \param order  Quadrature order.
   *  \param mu     Absissca
   *  \param wt     Weights
   */
  void generate_parameters(size_t order, vec_dbl &mu, vec_dbl &wt);

};

} // end namespace detran_angle

#endif /* GAUSSLEGENDRE_HH_ */

//---------------------------------------------------------------------------//
//              end of GaussLegendre.hh
//---------------------------------------------------------------------------//
