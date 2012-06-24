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
 *
 *
 */
class Uniform: public QuadratureMOC
{

public:

  /*!
   *  \brief Constructor
   *  \param dim                  Problem dimension (only 2 supported for now)
   *  \param num_azimuths_octant  Number of azimuths per octant
   *  \param num_space            Number of tracks per azimuth
   *  \param num_polar            Number of polar angles in half space
   *  \param polar                Polar quadrature string identifier
   */
  Uniform(int dim,
          int num_azimuths_octant,
          int num_space,
          int num_polar,
          std::string polar);

  /// SP Constructor.
  static SP<QuadratureMOC> Create(int dim,
                                  int num_azimuths_octant,
                                  int num_space,
                                  int num_polar,
                                  std::string polar)
  {
    SP_quadrature p;
    p = new Uniform(dim, num_azimuths_octant, num_space, num_polar, polar);
    return p;
  }

  ~Uniform(){};

private:


};

} // end namespace detran

#endif // UNIFORM_HH_ 

//---------------------------------------------------------------------------//
//              end of file Uniform.hh
//---------------------------------------------------------------------------//
