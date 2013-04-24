//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_1D.hh
 *  @author robertsj
 *  @date   Jun 9, 2012
 *  @brief  Equation_SC_1D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_SC_1D_HH_
#define detran_EQUATION_SC_1D_HH_

#include "Equation.hh"

namespace detran
{

/**
 *  @class Equation_SC_1D
 *  @brief Step characteristic discretization in one dimension.
 *
 *  In 1D, the step characteristic approximation is the same for
 *  discrete ordinates and MOC.  The outgoing flux is defined
 *  \f[
 *      \psi_{out} = A\psi_{in}  + B Q   \, ,
 *  \f]
 *  and the average segment flux is
 *  \f[
 *      \bar{\psi} = \frac{1}{l} \Big ( B \psi_{in} +  C Q \Big )  \, ,
 *  \f]
 *  where
 *  \f[
 *      A = e^{-\Sigma_t \tau} \, ,
 *  \f]
 *  \f[
 *      B = \frac{1}{\Sigma_t} ( 1- A ) \, ,
 *  \f]
 *  and
 *  \f[
 *      C = \frac{l}{\Sigma_t} \Big( 1- \frac{1-A}{\tau} \Big ) \, ,
 *  \f]
 *  where \f$ l \f$ is the step length and
 *  \f$ \tau = \Sigma_t l  \f$ is optical path length.
 *
 *  @sa Equation_SC_MOC
 */
class TRANSPORT_EXPORT Equation_SC_1D : public Equation<_1D>
{

public:

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
  Equation_SC_1D(SP_mesh mesh,
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

  /// Cosine
  double d_mu;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SC_1D.i.hh"

#endif /* detran_EQUATION_SC_1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_1D.hh
//---------------------------------------------------------------------------//
