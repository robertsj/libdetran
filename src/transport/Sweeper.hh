//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER_HH_
#define SWEEPER_HH_

// Other libtran headers
#include "Boundary.hh"
#include "Equation.hh"
#include "Material.hh"
#include "Mesh.hh"
#include "Quadrature.hh"
#include "State.hh"
#include "SweepSource.hh"

// Other libtran utils headers
#include "Definitions.hh"
#include "SP.hh"
#include "InputDB.hh"

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Sweeper
 * \brief Sweeper for discrete ordinates problems.
 *
 * The within-group transport equation is
 * \f[
 *      \mathbf{L}\psi = Q \, ,
 * \f]
 * where \f$ \mathbf{L} \f$ is the streaming and collision operator and
 * \f$ Q \f$ is a discrete representation of all source contributions.
 *
 * To invert the operator \f$ \mathbf{L} \f$, we "sweep" over the mesh for all
 * angles,
 * which gives us updated angular fluxes in each cell.  Actually, the
 * flux *moments* are updated, while the discrete angular flux is
 * optionally stored.
 *
 * Relevant input database entries:
 *   - store_angular_flux [int]
 *   - equation [string]
 *
 */
//---------------------------------------------------------------------------//
template <class D>
class Sweeper
{

public:

  typedef SP<Sweeper>                       SP_sweeper;
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  typedef typename
      Equation<D>::SP_equation              SP_equation;
  typedef Boundary<D>                       Boundary_T;
  typedef typename Boundary_T::SP_boundary  SP_boundary;
  typedef typename
      BoundaryTraits<D>::value_type         boundary_flux_type;
  //
  typedef typename
      SweepSource<D>::SP_sweepsource        SP_sweepsource;
  //
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   */
  Sweeper(SP_input input,
          SP_mesh mesh,
          SP_material material,
          SP_quadrature quadrature,
          SP_state state,
          SP_boundary boundary,
          SP_sweepsource sweepsource);

  /*!
   *  \brief SP Constructor.
   */
  static detran::SP<Sweeper<D> >
  Create(detran::SP<InputDB>                 input,
         detran::SP<detran::Mesh>            mesh,
         detran::SP<detran::Material>        material,
         detran::SP<detran::Quadrature>      quadrature,
         detran::SP<detran::State>           state,
         detran::SP<detran::Boundary<D> >    boundary,
         detran::SP<detran::SweepSource<D> > sweepsource)
  {
    SP_sweeper p;
    p = new Sweeper(input, mesh, material, quadrature,
                    state, boundary, sweepsource);
    return p;
  }

  /*!
   *  \brief Sweep over all angles and space.
   *
   *  Note, if the angular flux is to be updated,
   *  it is done directly to via State.  Having
   *  sweep take the flux as an explicit argument
   *  allows various input types (e.g. Krylov
   *  vectors) without having to go through State.
   *
   *  Note, here the source is limited to a moments-based
   *  one.  We may want a discrete source later; that can
   *  be added with a third argument of psi type.
   *
   */
  inline void sweep(moments_type &phi)
  {/* ... */}

  /// Setup the equations for the group
  void setup_group(int g)
  {
    d_g = g;
    d_equation->setup_group(g);
  }

  /// Allows the psi update to occur whenever needed
  void set_update_psi(bool v)
  {
    d_update_psi = v;
  }

  bool is_valid() const
  {
    return true;
  }

private:

  /// Input database
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// State vectors
  SP_state d_state;

  /// Equation
  SP_equation d_equation;

  /// Boundary
  SP_boundary d_boundary;

  /// Sweep source
  SP_sweepsource d_sweepsource;

  /// Current group
  int d_g;

  /// Update the angular flux?
  bool d_update_psi;

  /// Match incident/outgoing side with octant
  vec3_int d_face_index;

  /// Adjoint problem?
  bool d_adjoint;

  /// Allocate template-specific items.
  void setup(SP_material material)
  {}

  /// Mesh sweeper indices. \todo Allow adjoint.
  inline int index(int o, int dim, int ijk);

};

// Constructor
template <class D>
Sweeper<D>::Sweeper(SP_input input,
                    SP_mesh mesh,
                    SP_material material,
                    SP_quadrature quadrature,
                    SP_state state,
                    SP_boundary boundary,
                    SP_sweepsource sweepsource)
  : d_input(input)
  , d_mesh(mesh)
  , d_quadrature(quadrature)
  , d_state(state)
  , d_boundary(boundary)
  , d_sweepsource(sweepsource)
  , d_update_psi(false)
  , d_adjoint(false)
{
  Require(d_input);
  Require(d_mesh);
  Require(d_quadrature);
  Require(d_state);
  Require(d_boundary);
  Require(d_sweepsource);
  Require(material);
  setup(material);
  Ensure(d_equation);
  // Check whether we keep psi.
  if (d_input->check("store_angular_flux"))
  {
    d_update_psi = d_input->get<int>("store_angular_flux");
    std::cout << "update psi found " << d_update_psi << std::endl;
  }
  std::cout << "update psi found " << d_update_psi << std::endl;
}

// Mesh sweeper indices
template <class D>
inline int Sweeper<D>::index(int o, int dim, int ijk)
{
  if (dim == 1)
  {
    if ((o == 0 or o == 3 or o == 4 or o == 7))
      return ijk;
    else
      return d_mesh->number_cells_x() - ijk - 1;
  }
  if (dim == 2)
  {
    if ((o == 0 or o == 1 or o == 4 or o == 5))
      return ijk;
    else
      return d_mesh->number_cells_y() - ijk - 1;
  }
  if (dim == 3)
  {
    if ((o == 0 or o == 1 or o == 2 or o == 3))
      return ijk;
    else
      return d_mesh->number_cells_z() - ijk - 1;
  }
}

} // end namespace detran

// Template specializations
#include "Sweeper1D.t.hh"
#include "Sweeper2D.t.hh"
#include "Sweeper3D.t.hh"

#endif /* SWEEPER_HH_ */
