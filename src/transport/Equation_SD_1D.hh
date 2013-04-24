//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_SD_1D.hh
 *  @author robertsj
 *  @date   Jun 9, 2012
 *  @brief  Equation_SD_1D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_SD_1D_HH_
#define detran_EQUATION_SD_1D_HH_

#include "Equation.hh"

namespace detran
{

/**
 *  @class Equation_SD_1D
 *  @brief Step difference discretization in one dimension.
 *
 *  The step difference approximation defines
 *  @f[
 *      \psi_{i,n} = \left\{
 *        \begin{array}{l l}
 *          \psi_{i+1/2,n}      & \quad \text{if $\mu_n > 0$} \\
 *          \psi_{i-1/2,n}      & \quad \text{if $\mu_n < 0$} \\
 *        \end{array} \right\}
 *  @f]
 *
 *  This is a first order method, but is strictly positive.
 *  Moreover, it essentially defines a "flat flux" within
 *  the cell, and so it is consistent for DGM.
 *
 */
class TRANSPORT_EXPORT Equation_SD_1D : public Equation<_1D>
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
  Equation_SD_1D(SP_mesh mesh,
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

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Equation_SD_1D.i.hh"

#endif /* detran_EQUATION_SD_1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_SD_1D.hh
//---------------------------------------------------------------------------//
