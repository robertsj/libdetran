//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_2D.hh
 *  @author robertsj
 *  @date   Apr 10, 2012
 *  @brief  Equation_SC_2D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SC_2D_HH_
#define EQUATION_SC_2D_HH_

#include "Equation.hh"

namespace detran
{

/**
 *  @class Equation_SC_2D
 *  @brief Step characteristic discretization in two dimensions.
 *
 *  Reference: Lathrop, K. D. "Spatial differencing of the Transport
 *             Equation: Positivity vs. Accuracy", J. Comp. Phys. 4,
 *             475-498 (1969)
 *
 */
class TRANSPORT_EXPORT Equation_SC_2D : public Equation<_2D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Equation<_2D> >    SP_equation;
  typedef Equation<_2D>::SP_material              SP_material;
  typedef Equation<_2D>::SP_mesh                  SP_mesh;
  typedef Equation<_2D>::SP_quadrature            SP_quadrature;
  typedef Equation<_2D>::moments_type             moments_type;
  typedef Equation<_2D>::angular_flux_type        angular_flux_type;
  typedef Equation<_2D>::face_flux_type           face_flux_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  Equation_SC_2D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 const bool update_psi);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Solve for the cell-center and outgoing edge fluxes.
  inline void solve(const size_t i,
                    const size_t j,
                    const size_t k,
                    moments_type &source,
                    face_flux_type &psi_in,
                    face_flux_type &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /// Setup the equations for a group.
  void setup_group(const size_t g);

  /// Setup the equations for an octant.
  void setup_octant(const size_t octant);

  /// Setup the equations for an angle.
  void setup_angle(const size_t angle);


private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// X-directed coefficient, \f$ \Delta_x / |\mu| \f$.
  detran_utilities::vec_dbl d_alpha;

  /// Y-directed coefficient, \f$ \Delta_y / |\eta|  \f$.
  detran_utilities::vec_dbl d_beta;

  /**
   *  \brief Approximate exponential.
   *
   *  For SC, Denovo uses a 7th order truncate expansion for
   *  exp about zero for x <= 0.025.  For that value, the
   *  approximation yields an error of about 1e-15 (found in Maple
   *  using high precision).  This makes sense, since the error
   *  of the series is bounded by \f$ x^n / n! \f$ for \f$ x < 1 \f$
   *  and where \f$ n \f$ is the first truncated order.  Keeping this
   *  as the desired precision,
   *
   *  Here, we do the same but allow for the approximation's use
   *  over a intervals of width 0.1 from x = 0 to 1.0, which should
   *  bound most use cases.
   *
   *  \todo Performance study to see what effect this actually has!
   *
   */
  inline double exp_appx(double x);


};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SC_2D.i.hh"

#endif /* EQUATION_SC_2D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_2D.hh
//---------------------------------------------------------------------------//
