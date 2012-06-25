//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Collocated.hh
 * \brief  Collocated 
 * \author Jeremy Roberts
 * \date   Jun 25, 2012
 */
//---------------------------------------------------------------------------//

#ifndef COLLOCATED_HH_
#define COLLOCATED_HH_

#include "QuadratureMOC.hh"

namespace detran
{

/*!
 *  \class Collocated
 *  \brief A MOC quadrature with a finite set of origins along a side.
 */
class Collocated: public QuadratureMOC
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
  Collocated(int dim,
             int num_azimuths_octant,
             int num_space,
             int num_polar,
             std::string polar);

  /// SP Constructor.
  static SP<Quadrature> Create(int dim,
                               int num_azimuths_octant,
                               int num_space,
                               int num_polar,
                               std::string polar)
  {
    SP_quadrature p;
    p = new Collocated(dim, num_azimuths_octant, num_space, num_polar, polar);
    return p;
  }

  ~Collocated(){};

private:


};

} // end namespace detran

#endif // COLLOCATED_HH_ 

//---------------------------------------------------------------------------//
//              end of file Collocated.hh
//---------------------------------------------------------------------------//
