//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper2DMOC.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper2DMOC class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER2DMOC_HH_
#define SWEEPER2DMOC_HH_

// Detran
#include "BoundaryMOC.hh"
#include "Sweeper.hh"
#include "MeshMOC.hh"
#include "QuadratureMOC.hh"
#include "TrackDB.hh"
#include "Track.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 * \class Sweeper2DMOC
 * \brief Sweeper for 2D MOC problems.
 *
 * \todo An option would be to have a mesh tracked in sweep constructor
 */
//---------------------------------------------------------------------------//
template <class EQ>
class Sweeper2DMOC: public Sweeper<_2D>
{

public:

  typedef SP<Sweeper2DMOC>                  SP_sweeper;
  typedef Sweeper<_2D>                      Base;
  //
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef MeshMOC::SP_mesh                  SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef QuadratureMOC::SP_quadrature      SP_quadrature;
  //
  typedef EQ                                Equation_T;
  typedef BoundaryMOC<_2D>                  Boundary_T;
  typedef typename Boundary_T::SP_boundary  SP_boundary;
  //
  typedef typename
      SweepSource<_2D>::SP_sweepsource      SP_sweepsource;
  //
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;
  //
  typedef TrackDB::SP_trackdb               SP_trackdb;
  typedef Track::SP_track                   SP_track;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Tracked mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature for MOC.
   *  \param    state       State vectors.
   *  \param    boundary    Boundary based on tracks.
   *  \param    sweepsource Sweep source constructor.
   */
  Sweeper2DMOC(SP_input input,
               SP_mesh mesh,
               SP_material material,
               SP_quadrature quadrature,
               SP_state state,
               SP_boundary boundary,
               SP_sweepsource sweepsource)
  : Base(input,mesh,material,quadrature,
         state,sweepsource)
  , d_boundary(boundary)
  {
    d_tracks = mesh->tracks();
  }

  /// Virtual destructor
  virtual ~Sweeper2DMOC(){}

  /// SP Constructor
  static detran::SP<Sweeper2DMOC<EQ> >
  Create(detran::SP<InputDB>                    input,
         detran::SP<detran::Mesh>               mesh,
         detran::SP<detran::Material>           material,
         detran::SP<detran::Quadrature>         quadrature,
         detran::SP<detran::State>              state,
         detran::SP<detran::BoundaryBase<_2D> > boundary,
         detran::SP<detran::SweepSource<_2D> >  sweepsource)
  {
    SP_sweeper p(new Sweeper2DMOC(input, mesh, material, quadrature,
                                  state, boundary, sweepsource));
    return p;
  }

  /// Sweep.
  inline void sweep(moments_type &phi);

  /// Setup the equations for the group
  void setup_group(int g)
  {
    d_g = g;
  }

private:

  /// \name Private Data
  /// \{

  SP_boundary d_boundary;

  SP_trackdb d_tracks;

  /// \}

};

} // end namespace detran

#include "Sweeper2DMOC.i.hh"

#endif /* SWEEPER2DMOC_HH_ */
