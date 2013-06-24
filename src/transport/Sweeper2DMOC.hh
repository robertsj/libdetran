//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   Sweeper2DMOC.hh
 *  @author Jeremy Roberts
 *  @date   Mar 24, 2012
 *  @brief  Sweeper2DMOC class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_SWEEPER2DMOC_HH_
#define detran_SWEEPER2DMOC_HH_

#include "transport/Sweeper.hh"
#include "angle/ProductQuadrature.hh"
#include "boundary/BoundaryMOC.hh"
#include "geometry/Mesh.hh"
#include "geometry/TrackDB.hh"
#include "geometry/Track.hh"

namespace detran
{

/**
 *  @class Sweeper2DMOC
 *  @brief Sweeper for 2D MOC problems.
 */

template <class EQ>
class Sweeper2DMOC: public Sweeper<_2D>
{

public:
  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef detran_utilities::SP<Sweeper2DMOC>            SP_sweeper;
  typedef Sweeper<_2D>                                  Base;
  typedef typename Base::SP_state                       SP_state;
  typedef typename Base::SP_input                       SP_input;
  typedef typename Base::SP_material                    SP_material;
  typedef typename Base::SP_sweepsource                 SP_sweepsource;
  typedef typename Base::moments_type                   moments_type;
  typedef typename Base::angular_flux_type              angular_flux_type;
  typedef typename Base::SP_tally                       SP_tally;
  typedef typename Base::vec_int                        vec_int;
  typedef typename Base::vec2_int                       vec2_int;
  typedef typename Base::vec3_int                       vec3_int;
  typedef typename Base::size_t                         size_t;
  typedef EQ                                            Equation_T;
  typedef BoundaryMOC<_2D>                              Boundary_T;
  typedef typename Boundary_T::SP_boundary              SP_boundary;
  typedef detran_geometry::Mesh::SP_mesh                SP_mesh;
  typedef detran_angle::ProductQuadrature::SP_quadrature SP_quadrature;
  typedef detran_geometry::TrackDB::SP_trackdb          SP_trackdb;
  typedef detran_geometry::Track::SP_track              SP_track;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  /**
   *  @brief Constructor.
   *  @param    input       User input database.
   *  @param    mesh        Tracked mesh.
   *  @param    material    Material database.
   *  @param    quadrature  Angular quadrature for MOC.
   *  @param    state       State vectors.
   *  @param    boundary    Boundary based on tracks.
   *  @param    sweepsource Sweep source constructor.
   */
  Sweeper2DMOC(SP_input input,
               SP_mesh mesh,
               SP_material material,
               SP_quadrature quadrature,
               SP_state state,
               SP_boundary boundary,
               SP_sweepsource sweepsource);

  /// Virtual destructor
  virtual ~Sweeper2DMOC(){}

  /// SP Constructor
  static SP_sweeper
  Create(SP_input       input,
         SP_mesh        mesh,
         SP_material    material,
         SP_quadrature  quadrature,
         SP_state       state,
         SP_boundary    boundary,
         SP_sweepsource sweepsource);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL SWEEPERS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Sweep.
  inline void sweep(moments_type &phi);

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // MOC boundary
  SP_boundary d_boundary;
  // Track database
  SP_trackdb d_tracks;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE MEMBER DEFINITIONS
//---------------------------------------------------------------------------//

#include "Sweeper2DMOC.i.hh"

#endif /* detran_SWEEPER2DMOC_HH_ */
