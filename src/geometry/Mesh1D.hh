//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Mesh1D.hh
 * \author Jeremy Roberts
 * \date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef MESH1D_HH_
#define MESH1D_HH_

// Geometry headers
#include "Mesh.hh"

// System
#ifdef DETRAN_ENABLE_BOOST
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#endif

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Mesh1D
 *  \brief One-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class Mesh1D : public Mesh
{

public:

  typedef Mesh            Base;
  typedef Base::SP_mesh   SP_mesh;

  /*!
   *  \brief Constructor.
   *
   *  \param    xfm         Fine meshes per coarse mesh in x dimension.
   *  \param    xcme        Coarse mesh edges x dimension.

   *  \param    mat_map     Coarse mesh material map.
   */
  Mesh1D(vec_int xfm,  vec_dbl xcme,  vec_int mat_map);

  /*!
   *  \brief Constructor.
   *
   *  \param    xfme        Fine mesh edges x dimension.
   *  \param    mat_map     Fine mesh material map.
   */
   Mesh1D(vec_dbl xfme, vec_int mat_map);

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfm         Fine meshes per coarse mesh in x dimension.
    *  \param    xcme        Coarse mesh edges x dimension.
    *  \param    mat_map     Coarse mesh material map.
    */
   static SP<Mesh> Create(vec_int xfm,
                          vec_dbl xcme,
                          vec_int mat_map)
   {
     SP_mesh p;
     p = new Mesh1D(xfm, xcme, mat_map);
     return p;
   }

   /*!
    *  \brief SP Constructor.
    *
    *  \param    xfme        Fine mesh edges x dimension.
    *  \param    mat_map     Fine mesh material map.
    */
    static SP<Mesh> Create(vec_dbl xfme,  vec_int mat_map)
    {
      SP_mesh p;
      p = new Mesh1D(xfme, mat_map);
      return p;
    }

    bool is_valid() const
    {
      Ensure(d_number_cells_y == 1);
      Ensure(d_number_cells_z == 1);
      return true;
    }

protected:

  /*!
   *   We keep this as an option in the event inherited meshes need
   *   more flexibility.
   */
  Mesh1D() : Mesh(1) {}

private:

#ifdef DETRAN_ENABLE_BOOST
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)
  {
    ar & boost::serialization::base_object<Base>(*this);
  }
#endif

};

} // end namespace detran

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(detran::Mesh1D)
#endif

#endif /* MESH1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh1D.hh
//---------------------------------------------------------------------------//
