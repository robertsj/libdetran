//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   GaussSeidel.hh
 * \author robertsj
 * \date   Apr 9, 2012
 * \brief  GaussSeidel class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef GAUSSSEIDEL_HH_
#define GAUSSSEIDEL_HH_

// Detran
#include "MultigroupSolver.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class GaussSeidel
 * \brief Solves the multigroup transport equation via Gauss-Seidel.
 *
 * Relevant db entries:
 * - outer_norm_type (str) [default = "Linf"]
 */
//---------------------------------------------------------------------------//

template <class D>
class GaussSeidel: public MultigroupSolver<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<GaussSeidel<D> >     SP_solver;
  typedef MultigroupSolver<D>                       Base;
  typedef typename Base::SP_solver                  SP_base;
  typedef typename Base::SP_inner                   SP_inner;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_state                   SP_state;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_material                SP_material;
  typedef typename Base::SP_quadrature              SP_quadrature;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_externalsource          SP_externalsource;
  typedef typename Base::SP_fissionsource           SP_fissionsource;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *
   *  \param input             Input database.
   *  \param state             State vectors, etc.
   *  \param mesh              Problem mesh.
   *  \param mat               Material definitions.
   *  \param quadrature        Angular mesh.
   *  \param boundary          Boundary fluxes.
   *  \param external_source   User-defined external source.
   *  \param fission_source    Fission source.
   */
  GaussSeidel(SP_input           input,
              SP_state           state,
              SP_mesh            mesh,
              SP_material        material,
              SP_quadrature      quadrature,
              SP_boundary        boundary,
              SP_externalsource  q_e,
              SP_fissionsource   q_f);

  /// SP constructor
  static SP_solver
  Create(SP_input             input,
         SP_state             state,
         SP_mesh              mesh,
         SP_material          material,
         SP_quadrature        quadrature,
         SP_boundary          boundary,
         SP_externalsource    q_e,
         SP_fissionsource     q_f)
  {
    SP_solver p;
    p = new GaussSeidel(input, state, mesh, material,
                        quadrature, boundary, q_e, q_f);
    return p;
  }

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL MULTIGROUP SOLVERS MUST IMPLEMENT
  //-------------------------------------------------------------------------//

  /// Solve the multigroup equations.
  void solve();

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
  using Base::d_external_source;
  using Base::d_fissionsource;
  using Base::d_downscatter;
  using Base::d_number_groups;
  using Base::d_max_iters;
  using Base::d_tolerance;
  using Base::d_print_out;
  using Base::d_print_interval;
  using Base::d_inner_solver;

  /// Determines which norm to use (default is Linf)
  std::string d_norm_type;

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "GaussSeidel.i.hh"

#endif /* GAUSSSEIDEL_HH_ */
