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

// Utilities
#include "SP.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class SourceIteration
 * \brief 
 */
//===========================================================================//
template <class D>
class SourceIteration: public InnerIteration<D>
{

public:

  typedef detran_utils::SP<SourceIteration>     SP_inner;
  typedef InnerIteration<D>                     Base;
  typedef typename InnerIteration<D>::SP_inner  SP_base;
  // basic objects
  typedef detran_utils::InputDB::SP_input       SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename Boundary<D>::SP_boundary     SP_boundary;
  typedef typename MomentToDiscrete<D>::SP_MtoD SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;
  //
  typedef typename Sweeper<D>::SP_sweeper       SP_sweeper;
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

  /*!
   *  \brief SP Constructor.
   *
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    yfm         Fine meshes per coarse mesh in y dimension.
   *  \param    xcme        Coarse mesh edges x dimension.
   *  \param    ycme        Coarse mesh edges y dimension.
   *  \param    mat_map     Coarse mesh material map.
   */
  static detran_utils::SP<SourceIteration<D> >
  Create(detran_utils::SP<detran_utils::InputDB>   input,
         detran_utils::SP<detran::State>           state,
         detran_utils::SP<detran::Mesh>            mesh,
         detran_utils::SP<detran::Material>        material,
         detran_utils::SP<detran::Quadrature>      quadrature,
         detran_utils::SP<detran::Boundary<D> >    boundary,
         detran_utils::SP<detran::ExternalSource>  q_e,
         detran_utils::SP<detran::FissionSource>   q_f)
  {
    SP_inner p;
    p = new SourceIteration(input, state, mesh, material,
                            quadrature, boundary, q_e, q_f);
    return p;
  }

  /*!
   *  \brief Solve the within group equation.
   */
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
