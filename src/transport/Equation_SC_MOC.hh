//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SC_MOC.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_SC_MOC class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_SC_MOC_HH_
#define detran_EQUATION_SC_MOC_HH_

#include "Equation_MOC.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class Equation_SC_MOC
 *  @brief Step characteristic discretization for MOC
 *
 *  In the method of characteristics, the flux is solved for along a
 *  track assuming a flat source.  For a given incident flux into
 *  a track segment, we can define the outgoing segment flux
 *  @f[
 *      \psi_{out} = A\psi_{in}  + B Q   \, ,
 *  @f]
 *  and average segment flux
 *  @f[
 *      \bar{\psi} = \frac{1}{l} \Big ( B \psi_{in} +  C Q \Big )  \, ,
 *  @f]
 *  where
 *  @f[
 *      A = e^{-\Sigma_t \tau} \, ,
 *  @f]
 *  @f[
 *      B = \frac{1}{\Sigma_t} ( 1- A ) \, ,
 *  @f]
 *  and
 *  @f[
 *      C = \frac{l}{\Sigma_t} \Big( 1- \frac{1-A}{\tau} \Big ) \, ,
 *  @f]
 *  where \f$ l \f$ is the segment length and
 *  @f$ \tau = \Sigma_t l  \f$ is optical path length.
 *
 *  The step characteristic method is positive but only first-order
 *  accurate in space.
 *
 *  Reference:  A. Hebert, <em>Applied Reactor Physics</em>.
 *
 * \sa Equation_DD_MOC
 *
 */
//---------------------------------------------------------------------------//

class TRANSPORT_EXPORT Equation_SC_MOC: public Equation_MOC
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef Equation_MOC                      Base;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /// Constructor
  Equation_SC_MOC(SP_mesh mesh,
                  SP_material material,
                  SP_quadrature quadrature,
                  const bool update_psi);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MOC EQUATION TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Solve for the cell-center and outgoing edge fluxes.
  inline void solve(const size_t       region,
                    const double       length,
                    const double       width,
                    moments_type      &source,
                    double            &psi_in,
                    double            &psi_out,
                    moments_type      &phi,
                    angular_flux_type &psi);


  /// Setup the equations for a group.
  void setup_group(const size_t g);

  /// Setup the equations for an octant.
  void setup_octant(const size_t o);

  /// Setup the equations for an azimuth.
  void setup_azimuth(const size_t a);

  /// Setup the equations for a polar angle.
  void setup_polar(const size_t p);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Angular quadrature weights for an azimuth and all polar angles.
  detran_utilities::vec_dbl d_weights;

  /// Inverse polar sines
  detran_utilities::vec_dbl d_inv_sin;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SC_MOC.i.hh"

#endif /* detran_EQUATION_SC_MOC_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_MOC.hh
//---------------------------------------------------------------------------//
