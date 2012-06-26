//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SourceIteration.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  SourceIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef SOURCEITERATION_HH_
#define SOURCEITERATION_HH_

// Detran
#include "InnerIteration.hh"

#include <iostream>

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class SourceIteration
 * \brief 
 */
//---------------------------------------------------------------------------//

template <class D>
class SourceIteration: public InnerIteration<D>
{

public:

  typedef SP<SourceIteration>                   SP_inner;
  typedef InnerIteration<D>                     Base;
  typedef typename InnerIteration<D>::SP_inner  SP_base;
  // basic objects
  typedef InputDB::SP_input                     SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename BoundaryBase<D>::SP_boundary SP_boundary;
  typedef typename MomentToDiscrete<D>::SP_MtoD SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;
  //
  typedef typename Sweeper<D>::SP_sweeper       SP_sweeper;
  typedef typename
      SweepSource<D>::SP_sweepsource            SP_sweepsource;
  //
  typedef State::moments_type                   moments_type;

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
  SourceIteration(SP_input           input,
                  SP_state           state,
                  SP_mesh            mesh,
                  SP_material        material,
                  SP_quadrature      quadrature,
                  SP_boundary        boundary,
                  SP_externalsource  q_e,
                  SP_fissionsource   q_f);

  /// SP Constructor
  static SP<SourceIteration<D> >
  Create(SP<detran::InputDB>          input,
         SP<detran::State>            state,
         SP<detran::Mesh>             mesh,
         SP<detran::Material>         material,
         SP<detran::Quadrature>       quadrature,
         SP<detran::BoundaryBase<D> > boundary,
         SP<detran::ExternalSource>   q_e,
         SP<detran::FissionSource>    q_f)
  {
    SP_inner p(new SourceIteration(input, state, mesh, material,
                                   quadrature, boundary, q_e, q_f));
    return p;
  }

  /// Solve the within group equation.
  void solve(int g);

  // Make inherited data visible
  using Base::d_input;
  using Base::d_state;
  using Base::d_mesh;
  using Base::d_material;
  using Base::d_quadrature;
  using Base::d_boundary;
  using Base::d_sweeper;
  using Base::d_sweepsource;
  using Base::d_tolerance;
  using Base::d_max_iters;
  using Base::d_print_out;
  using Base::d_print_interval;
  //using Base::b_acceleration;

};

} // namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "SourceIteration.i.hh"

#endif /* SOURCEITERATION_HH_ */

//---------------------------------------------------------------------------//
//              end of SourceIteration.hh
//---------------------------------------------------------------------------//
