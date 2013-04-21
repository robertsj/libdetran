//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Mesh1D.hh
 *  @author Jeremy Roberts
 *  @date   Mar 23, 2012
 */
//---------------------------------------------------------------------------//

#ifndef detran_geometry_MESH1D_HH_
#define detran_geometry_MESH1D_HH_

#include "Mesh.hh"
#ifdef DETRAN_ENABLE_BOOST
#include <boost/serialization/access.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/export.hpp>
#endif

namespace detran_geometry
{

//---------------------------------------------------------------------------//
/**
 *  @class Mesh1D
 *  @brief One-dimensional Cartesian mesh.
 *
 *  This is mostly a convenience interface.
 */
//---------------------------------------------------------------------------//
class GEOMETRY_EXPORT Mesh1D : public Mesh
{

public:

  typedef Mesh Base;
  typedef Base::SP_mesh SP_mesh;

  /**
   *  @brief Constructor.
   *
   *  @param    xfm         Fine meshes per coarse mesh in x dimension.
   *  @param    xcme        Coarse mesh edges x dimension.
   *  @param    mat_map     Coarse mesh material map.
   */
  Mesh1D(vec_int xfm, vec_dbl xcme, vec_int mat_map);

  /**
   *  @brief Constructor.
   *
   *  @param    xfme        Fine mesh edges x dimension.
   *  @param    mat_map     Fine mesh material map.
   */
  Mesh1D(vec_dbl xfme, vec_int mat_map);

  /// SP constructor
  static detran_utilities::SP<Mesh>
  Create(vec_int xfm, vec_dbl xcme, vec_int mat_map)
  {
    SP_mesh p(new Mesh1D(xfm, xcme, mat_map));
    return p;
  }

  /// SP constructor
  static detran_utilities::SP<Mesh>
  Create(vec_dbl xfme, vec_int mat_map)
  {
    SP_mesh p(new Mesh1D(xfme, mat_map));
    return p;
  }

protected:

  /**
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

} // end namespace detran_geometry

#ifdef DETRAN_ENABLE_BOOST
BOOST_CLASS_EXPORT_KEY(detran_geometry::Mesh1D)
#endif

#endif /* detran_geometry_MESH1D_HH_ */

//---------------------------------------------------------------------------//
//              end of Mesh1D.hh
//---------------------------------------------------------------------------//
