//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryBase.hh
 * \brief  BoundaryBase 
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYBASE_HH_
#define BOUNDARYBASE_HH_

#include "discretization/DimensionTraits.hh"
#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/InputDB.hh"

namespace detran
{

/*!
 *  \page boundary Boundaries
 *
 *  One of the chief motivations for developing Detran was for it to
 *  support development of a response matrix code.
 *
 *  The various flavors of the response matrix method all employ the
 *  following idea: decompose a global mesh into several local meshes,
 *  and connect the various regions by approximate boundary conditions.
 *  These approximate boundary conditions are often defined in terms
 *  of "response functions", which say something about the flux
 *  production within a local mesh  (a flux "response") or the angular
 *  current leaving the partial mesh (a current "response") due to some
 *  unit source at a boundary or within the medium.
 *
 *  To support generation of current responses, a careful treatment of
 *  boundary conditions and boundary flux information is warranted. To
 *  that end, we suggest the following structure for boundary
 *  indexing.
 *
 *  In general, boundary fluxes will be stored in a manner consistent with
 *  the mesh.  For example, for a 2-D problem, the J-1 entry of the
 *  east boundary corresponds to north-east corner of the mesh.
 *
 *  To define a spatial-ordering scheme that is consistent for response
 *  generation, we take the natural ordering to be

 *    West:       [J – 1 .. 0], [0 .. K – 1]
 *    East:       [0 .. J – 1], [0 .. K – 1]
 *    South:      [K – 1 .. 0], [0 .. I – 1]
 *    North:      [0 .. K – 1], [0 .. I – 1]
 *    Bottom:     [I – 1 .. 0], [0 .. J – 1]
 *    Top:        [0 .. I – 1], [0 .. J – 1]
 *
 *  Using the spatial_index method, one can access the elements in this
 *  order.  For MOC, a different approach will be required.
 */


/*!
 *  \class BoundaryBase
 *  \brief Base class for boundary flux containers.
 *
 *  Ideally, boundaries should be extremely flexible in terms of accessing
 *  the elements in a desired order, which is needed for response generation.
 *  For fine mesh problems, indexing is relatively straightforward given a
 *  regular grid.  For MOC problems, however, the indexing is less
 *  straightforward.
 *
 *  Ultimately, I would like a relatively compact interface applicable to
 *  SN, MOC, and any future approximations.  My guess is this would entail
 *  some form of iterator, or perhaps use of a Visitor class.
 */

template <class D>
class BoundaryBase
{

public:

  //-------------------------------------------------------------------------//
  // ENUMERATIONS
  //-------------------------------------------------------------------------//

  /// Incident or outgoing flux
  enum inout
  {
      IN,
      OUT
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<BoundaryBase>        SP_boundary;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_utilities::size_t                  size_t;

  //-------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   */
  BoundaryBase(SP_input        input,
               SP_mesh         mesh)
    : d_input(input)
    , d_mesh(mesh)
    , d_has_reflective(false)
    , d_is_reflective(2*D::dimension, false)
    , d_boundary_flux_size(2*D::dimension, 0)
    , d_g(0)
  {
    Require(d_input);
    Require(d_input->check("number_groups"));
    d_number_groups = input->get<int>("number_groups");
    Require(d_number_groups > 0);
  }

  /// Virtual destructor
  virtual ~BoundaryBase(){}

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY TYPES MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  \param  g   Group of current solve
   */
  virtual void set(const size_t g) = 0;

  /*!
   *  \brief Update the boundaries for each sweep.
   *
   *  This updates all incident boundary fluxes
   *  using the current outgoing boundary fluxes
   *  in a group.  What happens in the update is
   *  a function of the underlying boundary
   *  condition.
   *
   *  \param  g   Group of current solve
   */
  virtual void update(const size_t g) = 0;

  /*!
   *  \brief Update the boundaries for a single angle.
   *
   *  This is an alternative update that only updates
   *  the incident boundary flux for a particular
   *  angle.  When called within a sweep, this allows
   *  the most recent boundary fluxes to be used,
   *  in effect producing a Gauss-Seidel iteration.
   *
   *  \note This cannot be used for Krylov solvers.
   *
   *  \param  g   Group of current solve
   *  \param  o   Octant being swept
   *  \param  a   Angle within octant being swept
   */
  virtual void update(const size_t g, const size_t o, const size_t a) = 0;

  /*!
   *  \brief Clear the boundary container for a group.
   *
   *  In some cases, a client might require homogeneous
   *  boundaries, perhaps after a fixed boundary has
   *  been used to construct a right hand side for a
   *  Krylov solve.
   *
   *  \param  g   Group of current solve
   */
  virtual void clear(const size_t g) = 0;

  /*
   *  \brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  virtual void set_incident(const size_t g, double *v) = 0;

  /*
   *  \brief Get the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  virtual void get_incident(const size_t g, double *v) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /// Clear all groups.
  void clear()
  {
    for (size_t g = 0; g < d_number_groups; ++g) clear(g);
  }

  /// Does the boundary have any reflective conditions?
  bool has_reflective() const
  {
    return d_has_reflective;
  }

  /// Is a side reflective?
  bool is_reflective(const size_t side) const
  {
    Require(side < d_is_reflective.size());
    return d_is_reflective[side];
  }

  /// Number of boundary flux values in a group on a side
  size_t boundary_flux_size(const size_t side) const
  {
    Require(side < d_boundary_flux_size.size());
    return d_boundary_flux_size[side];
  }

  /*!
   *  \brief Compute spatial index along a surface given actual index
   *
   *  All
   *
   *  \param d      Dimension (0, 1, or 2)
   *  \param s      Surface index (0..2*D)
   *  \param ijk    Actual index i, j, or k on the mesh
   */
  size_t spatial_index(const size_t dim, const size_t s, const size_t ijk) const
  {
    // Preconditions
    Require(dim < D::dimension);
    Require(s < D::dimension * 2);
    Require(ijk < d_mesh->number_cells(dim));

    size_t ijk_new = ijk;

    if (dim == 1)
    {
      if (s == Mesh::SOUTH)  return d_mesh->number_cells_x() - ijk - 1;
    }
    if (dim == 2)
    {
      if (s == Mesh::WEST)   return d_mesh->number_cells_y() - ijk - 1;
    }
    if (dim == 3)
    {
      if (s == Mesh::BOTTOM) return d_mesh->number_cells_z() - ijk - 1;
    }

    return ijk_new;
  }

protected:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  /// Input
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Do I have any reflective conditions?  (Krylov support)
  bool d_has_reflective;

  /// Vector of is it reflective? (Krylov support)
  std::vector<bool> d_is_reflective;

  /// Size of the boundary flux on a side in one group.
  detran_utilities::vec_int d_boundary_flux_size;

  /// Number of groups
  size_t d_number_groups;

  /// Current group being solved
  size_t d_g;

};

} // end namespace detran

#endif // BOUNDARYBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryBase.hh
//---------------------------------------------------------------------------//
