//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   BoundaryBase.hh
 *  @brief  BoundaryBase class definition
 *  @author Copyright (C) 2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_BOUNDARYBASE_HH_
#define detran_BOUNDARYBASE_HH_

#include "boundary/boundary_export.hh"
#include "transport/DimensionTraits.hh"
#include "geometry/Mesh.hh"
#include "utilities/DBC.hh"
#include "utilities/Definitions.hh"
#include "utilities/InputDB.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class BoundaryBase
 *  @brief Base class for boundary flux containers.
 */
//---------------------------------------------------------------------------//

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

  /// Get or set
  enum getset
  {
      GET,
      SET
  };

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<BoundaryBase>        SP_boundary;
  typedef detran_utilities::InputDB::SP_input       SP_input;
  typedef detran_geometry::Mesh                     Mesh;
  typedef detran_geometry::Mesh::SP_mesh            SP_mesh;
  typedef detran_utilities::size_t                  size_t;
  typedef D                                         D_T;

  //-------------------------------------------------------------------------//
  // CONSTRUCTORS & DESTRUCTORS
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param input
   *  @param mesh
   */
  BoundaryBase(SP_input        input,
               SP_mesh         mesh)
    : d_input(input)
    , d_mesh(mesh)
    , d_has_reflective(false)
    , d_is_reflective(2*D::dimension, false)
    , d_has_vacuum(true)
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

  /**
   *  @brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for a solve, such as any
   *  external boundary source.
   *
   *  @param  g   Group of current solve
   */
  virtual void set(const size_t g) = 0;

  /**
   *  @brief Update the boundaries for each sweep.
   *
   *  This updates all incident boundary fluxes using the current outgoing
   *  boundary fluxes in a group.  What happens in the update is a function
   *  of the underlying boundary condition.
   *
   *  @param  g   Group of current solve
   */
  virtual void update(const size_t g) = 0;

  /**
   *  @brief Update the boundaries for a single angle.
   *
   *  This is an alternative update that only updates the incident boundary
   *  flux for a particular angle.  When called within a sweep, this allows
   *  the most recent boundary fluxes to be used, yielding a better iteration.
   *
   *  @note This cannot be used for Krylov solvers.
   *
   *  @param  g   Group of current solve
   *  @param  o   Octant being swept
   *  @param  a   Angle within octant being swept
   */
  virtual void update(const size_t g,
                      const size_t o,
                      const size_t a) = 0;

  /**
   *  @brief Clear the boundary container for a group.
   *
   *  In some cases, a client might require homogeneous boundaries, perhaps
   *  after a fixed boundary has been used to construct a right hand side
   *  for a Krylov solve.
   *
   *  @param  g   Group of current solve
   */
  virtual void clear(const size_t g) = 0;


  /**
   *  @brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  @param g  Group of current sweep
   *  @param v  Pointer to data used in Krylov solve
   */
  virtual void psi(const      size_t g,
                   double    *v,
                   const int  inout,
                   const int  gs,
                   bool       onlyref = true) = 0;

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Clear all groups.
  void clear()
  {
    for (size_t g = 0; g < d_number_groups; ++g) clear(g);
  }

  /**
   *  @brief Clear any fixed boundary conditions.
   */
  virtual void clear_bc()
  {
    /* ... */
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

  /// Does the boundary have any vacuum conditions?
  bool has_vacuum() const
  {
    return d_has_vacuum;
  }

  /// Number of boundary flux values in a group on a side
  size_t boundary_flux_size(const size_t side) const
  {
    Require(side < d_boundary_flux_size.size());
    return d_boundary_flux_size[side];
  }

  /// Display boundray information and contents
  virtual void display(bool inout) const
  {
    /* ... */
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
  detran_utilities::vec_bool d_is_reflective;
  /// Do I have any vacuum conditions? (Krylov support)
  bool d_has_vacuum;
  /// Size of the boundary flux on a side in one group.
  detran_utilities::vec_int d_boundary_flux_size;
  /// Number of groups
  size_t d_number_groups;
  /// Current group being solved
  size_t d_g;

};

} // end namespace detran

#endif // detran_BOUNDARYBASE_HH_

//---------------------------------------------------------------------------//
//              end of file BoundaryBase.hh
//---------------------------------------------------------------------------//
