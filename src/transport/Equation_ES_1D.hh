//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_ES_1D.hh
 *  @brief  Equation_ES_1D
 *  @author Jeremy Roberts
 *  @date   Jan 7, 2013
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_ES_1D_HH_
#define detran_EQUATION_ES_1D_HH_

#include "Equation.hh"

namespace detran
{

/**
 *  @class Equation_ES_1D
 *  @brief Explicit slope discretization in one dimension.
 *
 *  The explicit slope discretization uses a weighted
 *  difference equation of the form
 *  @f[
 *      \psi_{n,j} = \frac{1 + \alpha_{n,j}}{2} \psi_{n,j+\frac{1}{2}}
 *                 + \frac{1 - \alpha_{n,j}}{2} \psi_{n,j-\frac{1}{2}}
 *                 - \frac{\alpha_{n,j}}{2\Sigma_{tj}} T_{n, j} \, ,
 *  @f]
 *  where
 *  @f[
 *      T_{n, j} = \frac{1}{2} \left ( \Sigma_{s,j} \hat{\psi}_{n,j}
 *               + \hat{Q}_{n, j} \right ) \, ,
 *  @f]
 *  contains slopes for the scattering and external source; that is,
 *  the full right hand side is a linear function.
 *
 *  The weights @f$ \alpha_{n,j} @f$ can be selected in several
 *  ways, including that of the SC method,
 *  @f[
 *      \alpha_{n,j} = \coth{ \tau_{n, j} } - \frac{1}{\tau_{n,j}} \,
 *  @f]
 *  and
 *  @f[
 *      \alpha_{n,j} = \frac{\tau_{n,j}}{\theta + |\tau_{n,j}|} \, ,
 *  @f]
 *  where
 *  @f[
 *      \tau_{n, j} = \frac{ \Sigma_{t, j} \Delta_j }{2 \mu_n} \, .
 *  @f]
 *  In the latter, @f$ \theta = 3 @f$ corresponds to a linear
 *  discontinuous, and @f$ \theta = 1 @f$ to lumped linear
 *  discontinuous.  Neither of these is implemented.
 *
 *  The flux slope is approximated as
 *
 *
 */
class Equation_ES_1D : public Equation<_1D>
{

public:

  /// Number of unknowns per cell
  static const int number_cell_unknowns = 1;

  /// Number of unknowns per surface
  static const int number_surface_unknowns = 1;

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Equation<_1D> >  SP_equation;
  typedef Equation<_1D>::SP_material            SP_material;
  typedef Equation<_1D>::SP_mesh                SP_mesh;
  typedef Equation<_1D>::SP_quadrature          SP_quadrature;
  typedef Equation<_1D>::moments_type           moments_type;
  typedef Equation<_1D>::angular_flux_type      angular_flux_type;
  typedef Equation<_1D>::face_flux_type         face_flux_type;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  Equation_ES_1D(SP_mesh mesh,
                 SP_material material,
                 SP_quadrature quadrature,
                 bool update_psi);

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

  /// X-directed coefficient, \f$ 2|\mu|/\Delta_x \f$.
  detran_utilities::vec_dbl d_coef_x;

};

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_ES_1D.i.hh"

} // end namespace detran

#endif // EQUATION_ES_1D_HH_ 

//---------------------------------------------------------------------------//
//              end of file Equation_ES_1D.hh
//---------------------------------------------------------------------------//
