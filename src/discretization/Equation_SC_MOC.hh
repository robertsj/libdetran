//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Equation_SC_MOC.hh
 * \author Jeremy Roberts
 * \date   Mar 31, 2012
 * \brief  Equation_SC_MOC class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef EQUATION_SC_MOC_HH_
#define EQUATION_SC_MOC_HH_

// Detran headers
#include "Equation_MOC.hh"

// Detran utilities
#include "Definitions.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Equation_SC_MOC
 *  \brief Step characteristic discretization for MOC
 *
 *  In the method of characteristics, the flux is solved for along a
 *  track assuming a flat source.  For a given incident flux into
 *  a track segment, we can define the outgoing segment flux
 *  \f[
 *      \psi_{out} = A\psi_{in}  + B Q   \, ,
 *  \f]
 *  and average segment flux
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
 *  where \f$ l \f$ is the segment length and
 *  \f$ \tau = \Sigma_t l  \f$ is optical path length.
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

class Equation_SC_MOC: public Equation_MOC
{

public:

  typedef Material::SP_material             SP_material;
  typedef MeshMOC::SP_mesh                  SP_mesh;
  typedef TrackDB::SP_trackdb               SP_trackdb;
  typedef QuadratureMOC::SP_quadrature      SP_quadrature;
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

  /*!
   *  \brief Constructor
   */
  Equation_SC_MOC(SP_mesh mesh,
                  SP_material material,
                  SP_quadrature quadrature,
                  bool update_psi);

  /// \name Public Interface
  /// \{

  /*!
   *   \brief Solve for the cell-center and outgoing edge fluxes.
   *
   *   See \ref Equation for full description.
   */
  inline void solve(int region,
                    double length,
                    moments_type &source,
                    double &psi_in,
                    double &psi_out,
                    moments_type &phi,
                    angular_flux_type &psi);


  /*!
   *  \brief Setup the equations for a group.
   *  \param g     Current group.
   */
  void setup_group(int g);

  /*!
   *  \brief Setup the equations for an octant.
   *  \param o    Current octant index.
   */
  void setup_octant(int o);

  /*!
   *  \brief Setup the equations for an azimuth.
   *  \param a    Azimuth within octant.
   */
  void setup_azimuth(int a);

  /*!
   *  \brief Setup the equations for a polar angle.
   *  \param p    Polar index.
   */
  void setup_polar(int p);

  /// \}

private:

  /// \name Private Data
  /// \{

  /// Angular quadrature weights for an azimuth and all polar angles.
  vec_dbl d_weights;

  /// Inverse polar sines
  vec_dbl d_inv_sin;

  /// \}

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SC_MOC.i.hh"

#endif /* EQUATION_SC_MOC_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SC_MOC.hh
//---------------------------------------------------------------------------//
