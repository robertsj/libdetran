//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   BoundaryMOC.hh
 * \brief  BoundaryMOC class definition
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef BOUNDARYMOC_HH_
#define BOUNDARYMOC_HH_

// Detran
#include "BoundaryBase.hh"
#include "BoundaryConditionMOC.hh"
#include "MeshMOC.hh"
#include "QuadratureMOC.hh"

// Utilities
#include "DBC.hh"

namespace detran
{

/*
 *  \class BoundaryMOC
 *  \brief Boundary flux container for MOC problems.
 *
 *  The method of characteristics solves the transport equation
 *  by sweeping along fixed tracks crossing the domain.  Tracks
 *  begin and end at a global boundary.  This class stores the
 *  angular flux for each track.
 *
 */

template <class D>
class BoundaryMOC : public BoundaryBase<D>
{

public:

  typedef SP<BoundaryMOC>               SP_boundary;
  typedef InputDB::SP_input             SP_input;
  typedef SP<MeshMOC>                   SP_mesh;
  typedef SP<QuadratureMOC>             SP_quadrature;
  typedef BoundaryConditionMOC<D>       BC_T;
  typedef typename BC_T::SP_bc          SP_bc;
  typedef BoundaryBase<D>               Base;
  typedef typename Base::SP_boundary    SP_base;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    quadrature  Angular quadrature.
   */
  BoundaryMOC(SP_input        input,
              SP_mesh         mesh,
              SP_quadrature   quadrature);

  /// SP Constructor.
  static SP<detran::BoundaryBase<D> >
  Create(SP<detran::InputDB>       input,
         SP<detran::Mesh>          mesh,
         SP<detran::Quadrature>    quadrature)
  {
    SP_boundary p(new BoundaryMOC(input, mesh, quadrature));
    return p;
  }

  /// \name Inherited Interface
  /// \{

  /*!
   *  \brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  \param  g   Group of current solve
   */
  void set(int g);

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
  void update(int g);

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
  void update(int g, int o, int a);

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
  void clear(int g);

  /*
   *  \brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is is support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void set_incident(int g, double *v)
  {
    THROW("IMPLEMENT ME");
  }

  /*
   *  \brief Get the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  void get_incident(int g, double *v)
  {
    THROW("IMPLEMENT ME");
  }

  /// \}

  /// DBC function
  bool is_valid() const
  {
    return true;
  }

private:

  // Expose base class members.
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_quadrature;
  using Base::d_number_groups;
  using Base::d_has_reflective;
  using Base::d_is_reflective;
  using Base::d_boundary_flux_size;
  using Base::inout;

  /// \name Private Data
  /// \{

  /// Boundary flux (energy, angle).(space^D)
  vec3_dbl d_boundary_flux;

  /// Vector of boundary conditions.
  std::vector<SP_bc> d_bc;

  /// \}

  /// \name Implementation
  /// \{

  /// Size the boundaries, etc.
  void initialize();

  /// \}

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "BoundaryMOC.i.hh"

#endif // BOUNDARYMOC_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryMOC.hh
//---------------------------------------------------------------------------//
