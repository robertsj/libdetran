//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FixedBoundary.hh
 *  @author robertsj
 *  @date   Jan 30, 2013
 *  @brief  FixedBoundary class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_FIXEDBOUNDARY_HH_
#define detran_FIXEDBOUNDARY_HH_

#include "BoundaryCondition.hh"

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class FixedBoundary
 *  @brief Fixed boundary condition for SN calculations.
 *
 *  This class stores a copy of the boundary container that
 *  can be directly filled by the user.
 */
//---------------------------------------------------------------------------//

template <class D>
class FixedBoundary : public BoundaryCondition<D>
{

public:

  //-------------------------------------------------------------------------//
  // TYPEDEFS
  //-------------------------------------------------------------------------//

  typedef BoundaryCondition<D>                      Base;
  typedef typename Base::SP_bc                      SP_bc;
  typedef typename Base::Boundary_T                 Boundary_T;
  typedef typename Base::SP_boundary                SP_boundary;
  typedef typename Base::SP_input                   SP_input;
  typedef typename Base::SP_mesh                    SP_mesh;
  typedef typename Base::SP_quadrature              SP_quadrature;
  typedef typename Base::size_t                     size_t;
  typedef typename BoundaryTraits<D>::value_type    bf_type;
  typedef typename std::vector<bf_type>             vec1_bflux;
  typedef typename std::vector<vec1_bflux>          vec2_bflux;
  typedef typename std::vector<vec2_bflux>          vec3_bflux;

  //-------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //-------------------------------------------------------------------------//

  FixedBoundary(SP_boundary boundary,
                const size_t side,
                SP_input input,
                SP_mesh mesh,
                SP_quadrature quadrature);

  //-------------------------------------------------------------------------//
  // ABSTRACT INTERFACE -- ALL BOUNDARY CONDITIONS MUST IMPLEMENT THESE
  //-------------------------------------------------------------------------//

  /// Set initial and/or fixed boundary condition.
  void set(const size_t g);

  /// Update a boundary following a sweep. FixedBoundary does nothing.
  void update(const size_t g){}

  /// Update a boundary for a given angle following a sweep.
  void update(const size_t g, const size_t o, const size_t a){}

  //-------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //-------------------------------------------------------------------------//

  /// Const access to a boundary flux using cardinal indices.
  const bf_type&
  operator()(const size_t o,
             const size_t a,
             const size_t g) const;

  /// Mutable access to boundary flux using cardinal indices.
  bf_type&
  operator()(const size_t o,
             const size_t a,
             const size_t g);

  /// Clear any currently set conditions
  void clear();

private:

  //-------------------------------------------------------------------------//
  // DATA
  //-------------------------------------------------------------------------//

  // Make inherited data visible
  using Base::d_boundary;
  using Base::d_side;
  using Base::d_input;
  using Base::d_mesh;
  using Base::d_quadrature;

  size_t d_number_groups;
  vec3_bflux d_psi;

};

} // end namespace detran

//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//

#include "FixedBoundary.i.hh"


#endif /* detran_FIXEDBOUNDARY_HH_ */
