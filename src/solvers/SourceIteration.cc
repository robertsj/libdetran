//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   SourceIteration.cc
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  SourceIteration class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

// Detran
#include "SourceIteration.hh"

namespace detran
{

template <class D>
SourceIteration<D>::SourceIteration(SP_input          input,
                                    SP_state          state,
                                    SP_mesh           mesh,
                                    SP_material       material,
                                    SP_quadrature     quadrature,
                                    SP_MtoD           MtoD,
                                    SP_boundary       boundary,
                                    SP_externalsource q_e,
                                    SP_fissionsource  q_f)
  :  InnerIteration<D>::InnerIteration(input, state, mesh, material, quadrature, MtoD, boundary, q_e, q_f)
{
  /* ... */
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of SourceIteration.cc
//---------------------------------------------------------------------------//
