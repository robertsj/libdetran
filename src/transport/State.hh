//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   State.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  State class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef STATE_HH_
#define STATE_HH_

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class State
 * \brief Problem state.
 *
 * The state of the problem at any point in the solution process can be
 * defined by the unknown flux moments (scalar flux and higher moments),
 * as these quantities are sufficient to describe reaction rates, which is
 * typically what we need (e.g. doses or fission rates).  For eigenvalue
 * problems, keff is also important.
 *
 * We also include the option to store the angular flux for each cell for
 * use in debugging and other applications.  This can be either the cell
 * centered angular or the cell faces angular flux, the latter of which
 * is required for CMFD.
 *
 * The moments are stored in \ref Moments_Field objects, with one such
 * object for each group stored in a vector.  The Moments_Field class
 * is templated on equation (for its number of cell unknowns) and
 * dimension (for dimension-specific indexing).
 *
 * The angular flux, when stored, is arranged similarly in a \ref
 * Psi_Cell_Field or Psi_Face_Field.  These are \em not yet
 * implemented.
 *
 * When needed, the underlying data array for a group can be accessed as
 * demonstrated by the following (I <b>think!</b>).
 * \code
 std::vector<Moments_Field<D> > group_moments;
 // build moments...
 double *localgrouparray;
 localgrouparray = &moments[g][0];
 * \endcode
 * This would be useful e.g. for filling a PETSc Vec object to use in one
 * of their Krylov schemes (rather than copying into another array).
 *
 * \todo Test whether referencing a Moments_Field object at [0] actually
 *       yields the array underneath.
 */
//---------------------------------------------------------------------------//

// Other libtran headers
#include "Definitions.hh"
#include "Quadrature.hh"


class State
{

  typedef detran_utils::vec2_dbl    moments_type;
  typedef detran_utils::vec3_dbl    angular_flux_type;
  typedef Mesh::SP_mesh             SP_mesh;
  //typedef Boundary::SP_boundary     SP_boundary;

public:

  State(SP_input        input,
        SP_mesh         mesh,
        SP_quadrature   quadrature);


private:

  /// Cell-center scalar flux moments, [space-moment, energy]
  moments_type d_moments;

  /// Cell-center angular flux, [space, angle, energy]
  angular_flux_type d_angular_flux;

  /// k-eigenvalue
  double d_eigenvalue;

  /// Mesh
  SP_mesh d_mesh;

  /// Boundary fluxes
  //SP_boundary d_boundary;

  /// Number of energy groups
  int d_number_groups;

};

} // end namespace detran

#endif /* STATE_HH_ */

//---------------------------------------------------------------------------//
//              end of State.hh
//---------------------------------------------------------------------------//
