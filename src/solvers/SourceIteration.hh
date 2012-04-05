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
  typedef InnerIteration<D>            Base;
  // basic objects
  typedef detran_utils::InputDB::SP_input       SP_input;
  typedef State::SP_state                       SP_state;
  typedef Mesh::SP_mesh                         SP_mesh;
  typedef Material::SP_material                 SP_material;
  typedef Quadrature::SP_quadrature             SP_quadrature;
  typedef typename Boundary<D>::SP_boundary     SP_boundary;
  typedef MomentToDiscrete::SP_MtoD             SP_MtoD;
  // source typedefs
  typedef ExternalSource::SP_source             SP_externalsource;
  typedef FissionSource::SP_source              SP_fissionsource;
  //
  typedef typename SweeperBase<D>::SP_sweeper   SP_sweeper;
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
  SourceIteration(SP_input          input,
                 SP_state           state,
                 SP_mesh            mesh,
                 SP_material        material,
                 SP_quadrature      quadrature,
                 SP_MtoD            MtoD,
                 SP_boundary        boundary,
                 SP_externalsource  q_e,
                 SP_fissionsource   q_f);

  /*!
   *  \brief Solve the within group equation.
   */
  virtual void solve(int g);

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
