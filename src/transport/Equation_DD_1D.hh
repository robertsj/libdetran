//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Equation_DD_1D.hh
 *  @author Jeremy Roberts
 *  @date   Mar 31, 2012
 *  @brief  Equation_DD_1D class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_EQUATION_DD_1D_HH_
#define detran_EQUATION_DD_1D_HH_

#include "transport/transport_export.hh"
#include "transport/Equation.hh"

namespace detran
{

/**
 *  @class Equation_DD_1D
 *  @brief Diamond difference discretization in one dimension.
 *
 *  See \ref Equation_DD_3D for a general description of the
 *  diamond difference approximation.
 */
class TRANSPORT_EXPORT Equation_DD_1D : public Equation<_1D>
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
  Equation_DD_1D(SP_mesh mesh,
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

#include "Equation_DD_1D.i.hh"

#endif /* detran_EQUATION_DD_1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Equation_DD_1D.hh
//---------------------------------------------------------------------------//
