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

namespace detran_angle
{

/*!
 *  \class Uniform
 *  \brief Uniformly-spaced azimuthal quadrature set.
 */
class ANGLE_EXPORT Uniform: public QuadratureMOC
{

public:

  //-------------------------------------------------------------------------//
  // PUBLIC INTERFACE
  //-------------------------------------------------------------------------//

  /*!
   *  \brief Constructor
   *  \param dim                  Problem dimension (only 2 supported for now)
   *  \param num_azimuths_octant  Number of azimuths per octant
   *  \param num_space            Number of tracks per azimuth
   *  \param num_polar            Number of polar angles in half space
   *  \param polar                Polar quadrature string identifier
   */
  Uniform(size_t dim,
          size_t num_azimuths_octant,
          size_t num_space,
          size_t num_polar,
          std::string polar);

  /// SP Constructor.
  static detran_utilities::SP<Quadrature>
    Create(size_t dim,
           size_t num_azimuths_octant,
           size_t num_space,
           size_t num_polar,
           std::string polar)
  {
    SP_quadrature
      p(new Uniform(dim, num_azimuths_octant, num_space, num_polar, polar));
    return p;
  }

  ~Uniform(){};

private:

};

} // end namespace detran_angle

#endif // UNIFORM_HH_ 

//---------------------------------------------------------------------------//
//              end of file Uniform.hh
//---------------------------------------------------------------------------//
