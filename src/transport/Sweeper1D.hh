//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Sweeper1D.hh
 * \author Jeremy Roberts
 * \date   Mar 24, 2012
 * \brief  Sweeper1D class definition.
 * \note   Copyright (C) 2012 Jeremy Roberts.
 */
//---------------------------------------------------------------------------//

#ifndef SWEEPER1D_HH_
#define SWEEPER1D_HH_

// Detran
#include "Boundary.hh"
#include "Sweeper.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/*!
 *  \class Sweeper1D
 *  \brief Sweeper for 1D discrete ordinates problems.
 */
//---------------------------------------------------------------------------//
template <class EQ>
class Sweeper1D: public Sweeper<_1D>
{

public:

  typedef SP<Sweeper1D>                     SP_sweeper;
  typedef Sweeper<_1D>                      Base;
  //
  typedef State::SP_state                   SP_state;
  typedef InputDB::SP_input                 SP_input;
  typedef Mesh::SP_mesh                     SP_mesh;
  typedef Material::SP_material             SP_material;
  typedef Quadrature::SP_quadrature         SP_quadrature;
  //
  typedef EQ                                Equation_T;
  typedef Boundary<_1D>                     Boundary_T;
  typedef typename Boundary_T::SP_boundary  SP_boundary;
  typedef typename
      BoundaryTraits<_1D>::value_type       boundary_flux_type;
  //
  typedef typename
      SweepSource<_1D>::SP_sweepsource      SP_sweepsource;
  //
  typedef State::moments_type               moments_type;
  typedef State::angular_flux_type          angular_flux_type;

  /*!
   *  \brief Constructor.
   *
   *  \param    input       User input database.
   *  \param    mesh        Cartesian mesh.
   *  \param    material    Material database.
   *  \param    quadrature  Angular quadrature.
   *  \param    state       State vectors.
   *  \param    boundary    Boundary based on mesh.
   *  \param    sweepsource Sweep source constructor.
   */
  Sweeper1D(SP_input input,
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
    Require(boundary);
  }

  /// Virtual destructor
  virtual ~Sweeper1D(){}

  /// SP Constructor
  static detran::SP<Sweeper1D<EQ> >
  Create(detran::SP<InputDB>                    input,
         detran::SP<detran::Mesh>               mesh,
         detran::SP<detran::Material>           material,
         detran::SP<detran::Quadrature>         quadrature,
         detran::SP<detran::State>              state,
         detran::SP<detran::BoundaryBase<_1D> > boundary,
         detran::SP<detran::SweepSource<_1D> >  sweepsource)
  {
    SP_sweeper p(new Sweeper1D(input, mesh, material, quadrature,
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

  SP_boundary d_boundary;

};

} // end namespace detran

#include "Sweeper1D.i.hh"

#endif /* SWEEPER1D_HH_ */
