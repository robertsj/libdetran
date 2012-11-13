//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   GaussSeidelMG.hh
 *  @author robertsj
 *  @date   Apr 9, 2012
 *  @brief  GaussSeidelMG class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_MGSOLVERGS_HH_
#define detran_MGSOLVERGS_HH_

#include "MGTransportSolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class GaussSeidelMG
 *  @brief Solves the multigroup transport equation via Gauss-Seidel.
 *
 *  Relevant db entries:
 *  - outer_norm_type (str) [default = "Linf"]
 */
//---------------------------------------------------------------------------//

template <class D>
class MGSolverGS: public MGTransportSolver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef MGTransportSolver<D>                      Base;
  typedef typename Base::SP_solver                  SP_solver;
  typedef typename Base::SP_wg_solver               SP_wg_solver;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_quadrature              SP_quadrature;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::vec_externalsource         vec_externalsource;
  typedef typename Base::SP_fissionsource           SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor
   *  @param state             State vectors, etc.
   *  @param material          Material definitions.
   *  @param boundary          Boundary fluxes.
   *  @param q_e               Vector of user-defined external sources
   *  @param q_f               Fission source.
   *  @param multiply          Flag for a multiplying fixed source problem
   */
  MGSolverGS(SP_state                   state,
             SP_material                material,
             SP_boundary                boundary,
             const vec_externalsource  &q_e,
             SP_fissionsource           q_f,
             bool                       multiply = false);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  void solve(const double keff = 1.0);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Expose base members.
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_quadrature;
  using Base::d_boundary;
  using Base::d_externalsources;
  using Base::d_fissionsource;
  using Base::d_downscatter;
  using Base::d_number_groups;
  using Base::d_maximum_iterations;
  using Base::d_tolerance;
  using Base::d_print_level;
  using Base::d_print_interval;
  using Base::d_adjoint;
  using Base::d_wg_solver;
  using Base::d_multiply;

  /// Determines which norm to use (default is Linf)
  std::string d_norm_type;

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "MGSolverGS.i.hh"

#endif /* detran_MGSOLVERGS_HH_ */
