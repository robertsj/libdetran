//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Boundary.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Boundary class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARY_HH_
#define BOUNDARY_HH_

// Detran
#include "Quadrature.hh"
#include "Mesh.hh"
#include "Material.hh"
#include "Traits.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "InputDB.hh"
#include "SP.hh"

namespace detran
{

/*!
 *  \brief Boundary traits for simplify type access.
 */
template <class D>
class BoundaryTraits
{
public:
  typedef detran_utils::vec2_dbl value_type;
};
template <>
class BoundaryTraits<_2D>
{
public:
  typedef detran_utils::vec_dbl value_type;
};
template <>
class BoundaryTraits<_1D>
{
public:
  typedef double value_type;
};

//---------------------------------------------------------------------------//
/*!
 * \class Boundary
 * \brief Boundary flux container.
 *
 * Since arbitrary boundary functions are integral to the response matrix
 * method, it helps to have a really easy interface for handling boundary
 * information.  The Boundary class contains all boundary fluxes (incident
 * and outgoing) for each surface.  The exact type of a boundary flux
 * is templated. E.g. a 1d boundary flux (for an angle and group) is
 * just a double, whereas for a 3d problem, it is a 2-d array.
 *
 * Keeping the fundamental type templated should allow this to be used in
 * all the SN stuff and potentially MOC applications. 
 *
 * All boundary fluxes are stored as follows:
 * for g = 0, ng
 *   for a = 0, na
 *     for i,j,k
 *
 * Note, the angles are order by the quadrature as follows:
 * for octant = 0, no
 *   for angle = 0, na in octant
 *
 * For most applications, the order of the mu/eta/ksi within an octant
 * doesn't matter.  However, for response generation, it helps a lot
 * to have a standard ordering.  Hence, we enforce
 * for all mu
 *  for all eta
 *    for all xi
 *
 */
//---------------------------------------------------------------------------//
template <class D>
class Boundary : public detran_utils::Object
{

public:

  enum inout
  {
      IN,
      OUT
  };

  typedef detran_utils::SP<Boundary>                SP_boundary;
  typedef detran_utils::InputDB::SP_input           SP_input;
  typedef Mesh::SP_mesh                             SP_mesh;
  typedef Quadrature::SP_quadrature                 SP_quadrature;
  typedef typename BoundaryTraits<D>::value_type    boundary_flux_type;
  typedef std::vector<boundary_flux_type>           vec_boundary_flux;
  typedef std::vector<vec_boundary_flux>            vec2_boundary_flux;
  typedef std::vector<vec2_boundary_flux>           vec3_boundary_flux;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    quadrature  Angular quadrature.
   */
  Boundary(SP_input        input,
           SP_mesh         mesh,
           SP_quadrature   quadrature);

  /*!
   *  \brief SP Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    quadrature  Angular quadrature.
   */
  static SP_boundary Create(detran_utils::SP<detran_utils::InputDB> input,
                            detran_utils::SP<detran::Mesh>          mesh,
                            detran_utils::SP<detran::Quadrature>    quadrature)
  {
    SP_boundary p;
    p = new Boundary(input, mesh, quadrature);
    return p;
  }

  /*
   *  \brief Const access to a boundary flux using cardinal indices.
   *
   *  This (and the mutable version) interface is for use
   *  in sweeping, where octants and angles are cycled. 
   *
   *  \param    side  Side index.
   *  \param    o     Octant index.
   *  \param    a     Angle index (within octant).
   *  \param    g     Energy group.
   *  \return         Constant reference to boundary flux.
   */
  const boundary_flux_type& 
  operator()(int side, int o, int a, int g) const;

  /*
   *  \brief Mutable access to boundary flux using cardinal indices.
   */
  boundary_flux_type& operator()(int side, int o, int a, int g);

  /*
   *  \brief Const access to ordered incident flux.
   */
  const boundary_flux_type& incident(int side, int angle, int g) const;

  /*
   *  \brief Mutable access to incident boundary flux.
   */
  boundary_flux_type& incident(int side, int angle, int g);

  /*
   *  \brief Const access to ordered outgoing flux.
   */
  const boundary_flux_type& outgoing(int side, int angle, int g) const;

  /*
   *  \brief Mutable access to outgoing boundary flux.
   */
  boundary_flux_type& outgoing(int side, int angle, int g);

  /*!
   *  \brief  Map the local index to cardinal index.
   *
   *  In some cases, we need the boundary in its local order,
   *  e.g. left to right in space and angle.
   *
   *  \param    angle
   *
   */
  virtual int ordered_angle(int side, int angle, int inout) const;

  /// Return the input.
  SP_input get_input()
  {
    return d_input;
  }

  /// Return the mesh.
  SP_mesh get_mesh()
  {
    return d_mesh;
  }

  /// Return the quadratur.
  SP_quadrature get_quadrature()
  {
    return d_quadrature;
  }

  virtual bool is_valid() const
  {
    /* ... */
  }

private:

  /// \name Data
  /// \{

  /// Input database
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// Number of energy groups
  int d_number_groups;

  /// Boundary flux (side, energy, angle, space)
  vec3_boundary_flux d_boundary_flux;

  /// \}

  /// \name Implementation
  /// \{

  void initialize();

  /// \}

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "Boundary.i.hh"

#endif /* BOUNDARY_HH_ */

//---------------------------------------------------------------------------//
//              end of Boundary.hh
//---------------------------------------------------------------------------//
