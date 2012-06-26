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

// Detran
#include "Quadrature.hh"
#include "Mesh.hh"
#include "Traits.hh"

// Utilities
#include "DBC.hh"
#include "InputDB.hh"

namespace detran
{

/*!
 *  \class BoundaryBase
 *  \brief Base class for boundary flux containers.
 */

template <class D>
class BoundaryBase: public Object
{

public:

  /// Incident or outgoing flux index
  enum inout
  {
      IN,
      OUT
  };

  /// \name Useful Typedefs
  /// \{

  typedef SP<BoundaryBase>                          SP_boundary;
  typedef InputDB::SP_input                         SP_input;
  typedef Mesh::SP_mesh                             SP_mesh;
  typedef Quadrature::SP_quadrature                 SP_quadrature;

  /// \}

  /*!
   *  \brief Constructor
   */
  BoundaryBase(SP_input        input,
               SP_mesh         mesh,
               SP_quadrature   quadrature)
    : d_input(input)
    , d_mesh(mesh)
    , d_quadrature(quadrature)
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

  /// \name Public Interface
  /// \{

  /*!
   *  \brief Set the boundaries for a within-group solve.
   *
   *  This sets any boundaries that must be fixed for
   *  a solve, such as any external boundary source.
   *
   *  \param  g   Group of current solve
   */
  virtual void set(int g) = 0;

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
  virtual void update(int g) = 0;

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
  virtual void update(int g, int o, int a) = 0;

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
  virtual void clear(int g) = 0;

  /*
   *  \brief Set the entire group boundary flux for reflecting sides.
   *
   *  This is is support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  virtual void set_incident(int g, double *v) = 0;

  /*
   *  \brief Get the entire group boundary flux for reflecting sides.
   *
   *  This is in support of Krylov solvers.
   *
   *  \param g  Group of current sweep
   *  \param v  Pointer to data used in Krylov solve
   */
  virtual void get_incident(int g, double *v) = 0;

  /// \}

  /// Does the boundary have any reflective conditions?
  bool has_reflective() const
  {
    return d_has_reflective;
  }

  /// Is a side reflective?
  bool is_reflective(int side) const
  {
    Require(side < d_is_reflective.size());
    return d_is_reflective[side];
  }


protected:

  /// \name Protected Data
  /// \{

  /// Input
  SP_input d_input;

  /// Mesh
  SP_mesh d_mesh;

  /// Angular quadrature
  SP_quadrature d_quadrature;

  /// Do I have any reflective conditions?  (Krylov support)
  bool d_has_reflective;

  /// Vector of is it reflective? (Krylov support)
  std::vector<bool> d_is_reflective;

  /// Size of the boundary flux on a side in one group.
  vec_int d_boundary_flux_size;

  /// Number of groups
  u_int d_number_groups;

  /// Current group being solved
  u_int d_g;

  /// \}

private:

  /// \name Private Data
  /// \{




  /// \}

};


} // end namespace detran

#endif // BOUNDARYBASE_HH_ 

//---------------------------------------------------------------------------//
//              end of file BoundaryBase.hh
//---------------------------------------------------------------------------//
