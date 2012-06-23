//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Uniform.hh
 * \brief  Uniform class definition.
 * \author Jeremy Roberts
 * \date   Jun 22, 2012
 */
//---------------------------------------------------------------------------//

#ifndef UNIFORM_HH_
#define UNIFORM_HH_

#include "QuadratureMOC.hh"

namespace detran
{

/*!
 *  \class Uniform
 *  \brief Uniformly-spaced azimuthal quadrature set.
 *
 *  We enforce an octant-symmetric definition.
 *
 */
class Uniform: public QuadratureMOC
{

public:

  Uniform(int dim,
          int num_azimuths_octant,
          int num_space,
          int num_polar,
          std::string polar);

  ~Uniform(){};

private:


};

} // end namespace detran

#endif // UNIFORM_HH_ 

//---------------------------------------------------------------------------//
//              end of file Uniform.hh
//---------------------------------------------------------------------------//
