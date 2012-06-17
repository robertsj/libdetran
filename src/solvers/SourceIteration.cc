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
//
//#include "CMR.hh"

namespace detran
{

// Constructor
template <class D>
SourceIteration<D>::SourceIteration(SP_input          input,
                                    SP_state          state,
                                    SP_mesh           mesh,
                                    SP_material       material,
                                    SP_quadrature     quadrature,
                                    SP_boundary       boundary,
                                    SP_externalsource q_e,
                                    SP_fissionsource  q_f)
  :  InnerIteration<D>::InnerIteration(input,
                                       state,
                                       mesh,
                                       material,
                                       quadrature,
                                       boundary,
                                       q_e,
                                       q_f)
{
  //b_acceleration = new CMR<D>(mesh, material, quadrature);
  d_sweeper->set_update_boundary(true);
}

} // end namespace detran

//---------------------------------------------------------------------------//
//              end of SourceIteration.cc
//---------------------------------------------------------------------------//
