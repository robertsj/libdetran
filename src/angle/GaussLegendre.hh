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

  typedef detran_utils::SP<GaussLegendre> SP_quadrature;
  typedef detran_utils::SP<Quadrature>    SP_base;
  typedef Quadrature::vec_int             vec_int;
  typedef Quadrature::vec_dbl             vec_dbl;
  typedef Quadrature::vec2_int            vec2_int;
  typedef Quadrature::vec2_dbl            vec2_dbl;

  /*!
   *  \brief Constructor.
   *
   *  \param    order       Quadrature order.
   */
  GaussLegendre(int order);

  void display() const;

private:

  /*!
   *  \brief Compute GL parameters for arbitrary order.
   *  \param order  Quadrature order.
   *  \param mu     Absissca
   *  \param wt     Weights
   */
  void generate_parameters(int order, double *mu, double *wt);

};

} // end namespace detran

#endif /* GAUSSLEGENDRE_HH_ */

//---------------------------------------------------------------------------//
//              end of GaussLegendre.hh
//---------------------------------------------------------------------------//
