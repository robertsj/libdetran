//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   CurrentTally.hh
 * \brief  CurrentTally 
 * \author Jeremy Roberts
 * \date   Aug 7, 2012
 */
//---------------------------------------------------------------------------//

#ifndef CURRENTTALLY_HH_
#define CURRENTTALLY_HH_

// Detran
#include "CoarseMesh.hh"
#include "Equation.hh"
#include "Quadrature.hh"
#include "Traits.hh"

// Utilities
#include "Definitions.hh"
#include "SP.hh"

// System
#include <vector>

namespace detran
{

/*!
 *  \class CurrentTally
 *  \brief Records current or other angular data at cell boundaries
 *
 *
 */
template <class D>
class CurrentTally
{

public:

  typedef SP<CurrentTally>                  SP_currenttally;
  typedef CoarseMesh::SP_coarsemesh         SP_coarsemesh;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef CoarseMesh::SP_mesh               SP_mesh;
  typedef typename
    EquationTraits<D>::face_flux_type       face_flux_type;

  /*!
   *  \brief Constructor
   *  \param mesh
   *  \param quadrature
   *  \param number_groups
   */
  CurrentTally(SP_coarsemesh mesh,
               SP_quadrature quadrature,
               const u_int number_groups);

  /*!
   *  \brief Add angular flux to the current tally
   *
   *  \param  i           x mesh index
   *  \param  j           y mesh index
   *  \param  k           z mesh index
   *  \param  g           group index
   *  \param  o           octant
   *  \param  a           angle within octant
   *  \param  psi         edge angular flux
   *  \param  incident    flag indicating incident boundary surface
   */
  inline void tally(const u_int i,
                    const u_int j,
                    const u_int k,
                    const u_int g,
                    const u_int o,
                    const u_int a,
                    const face_flux_type psi,
                    const bool incident = false);

  /*!
   *  \brief Get the partial current from a surface and sense
   *  \param  i       x mesh index
   *  \param  j       y mesh index
   *  \param  k       z mesh index
   *  \param  g       group index
   *  \param  axis    0, 1, or 2 for x, y, or z
   *  \param  sense   true for positive (e.g. +x)
   */
  inline double partial_current(const u_int i,
                                const u_int j,
                                const u_int k,
                                const u_int g,
                                const u_int axis,
                                const u_int sense);

private:

  /// Coarse mesh
  SP_coarsemesh d_coarsemesh;

  /// Quadrature
  SP_quadrature d_quadrature;

  /// Number of groups
  const u_int d_number_groups;

  /// Partial currents [dimension][group][sense][index]
  std::vector<vec3_dbl> d_partial_current;

  /// Fine mesh coarse edge flags
  vec2_int d_coarse_edge_flag;

  /// Shifts the spatial index depending on direction
  vec2_int d_octant_shift;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "CurrentTally.i.hh"

#endif // CURRENTTALLY_HH_ 

//---------------------------------------------------------------------------//
//              end of file CurrentTally.hh
//---------------------------------------------------------------------------//
