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
 */
//---------------------------------------------------------------------------//

template <class D>
class GaussSeidel: public MultigroupSolver<D>
{

public:

  typedef SP<GaussSeidel<D> >                   SP_solver;
  typedef MultigroupSolver<D>                   Base;
  typedef typename Base::SP_solver              SP_base;
  typedef typename InnerIteration<D>::SP_inner  SP_inner;
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;

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

  /*!
   *  \brief SP Constructor.
   */
  static SP<GaussSeidel<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::ExternalSource>   q_e,
         SP<detran::FissionSource>    q_f)
  {
    SP_solver p;
    p = new GaussSeidel(input, state, mesh, material,
                        quadrature, boundary, q_e, q_f);
    return p;
  }

  /// Solve the multigroup equations.
  void solve();

  /// Unimplemented DBC function.
  bool is_valid() const
  {return true;}


private:

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

  /// \name Private Data
  /// \{



  /// \}

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "GaussSeidel.i.hh"

#endif /* GAUSSSEIDEL_HH_ */
