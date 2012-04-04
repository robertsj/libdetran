//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   FissionSource.hh
 * \author robertsj
 * \date   Apr 4, 2012
 * \brief  FissionSource class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts. 
 */
//---------------------------------------------------------------------------//

#ifndef FISSIONSOURCE_HH_
#define FISSIONSOURCE_HH_

// Detran
#include "Material.hh"
#include "Mesh.hh"

// Utilities
#include "DBC.hh"
#include "Definitions.hh"
#include "SP.hh"

namespace detran
{

//===========================================================================//
/*!
 * \class FissionSource
 * \brief 
 */
//===========================================================================//

class FissionSource
{

public:

  typedef detran_utils::SP<FissionSource>     SP_source;

};

} // end namespace detran

#endif /* FISSIONSOURCE_HH_ */

//---------------------------------------------------------------------------//
//              end of FissionSource.hh
//---------------------------------------------------------------------------//
