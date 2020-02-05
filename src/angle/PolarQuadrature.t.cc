//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PolarQuadrature.t.cc
 *  @brief PolarQuadrature explicit instantiations
 *  @note  Copyright (C) 2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "angle/PolarQuadrature.t.hh"

namespace detran_angle
{

template class PolarQuadrature<BaseGL>;
template class PolarQuadrature<BaseDGL>;
template class PolarQuadrature<BaseGC>;
template class PolarQuadrature<BaseDGC>;
template class PolarQuadrature<BaseUniform>;
template class PolarQuadrature<BaseUniformCosine>;
template class PolarQuadrature<BaseSimpson>;
template class PolarQuadrature<TabuchiYamamoto>;
template class PolarQuadrature<AbuShumaysDoubleRange>;
template class PolarQuadrature<TriGauss>;

} // end namespace detran_angle

//----------------------------------------------------------------------------//
//              end of PolarQuadrature.t.cc
//----------------------------------------------------------------------------//

