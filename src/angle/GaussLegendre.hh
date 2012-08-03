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

namespace detran
{

class GaussLegendre : public Quadrature
{

public:

  typedef SP<GaussLegendre> SP_quadrature;
  typedef SP<Quadrature>    SP_base;

  /*!
   *  \brief Constructor.
   *
   *  \param    order       Quadrature order.
   */
  GaussLegendre(int order);

  /*!
   *  \brief SP Constructor.
   *
   *  \param    order       Quadrature order.
   */
  static SP<Quadrature> Create(int order)
  {
    SP_quadrature p;
    p = new GaussLegendre(order);
    return p;
  }

  void display() const;

private:


  /// \name Implementation
  /// \{

  /*!
   *  \brief Compute GL parameters for arbitrary order.
   *  \param order  Quadrature order.
   *  \param mu     Absissca
   *  \param wt     Weights
   */
  void generate_parameters(int order, vec_dbl &mu, vec_dbl &wt);

  ///\}
};

} // end namespace detran

#endif /* GAUSSLEGENDRE_HH_ */

//---------------------------------------------------------------------------//
//              end of GaussLegendre.hh
//---------------------------------------------------------------------------//
