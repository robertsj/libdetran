//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   FixedBoundary.i.hh
 *  @author robertsj
 *  @date   Jan 30, 2013
 *  @brief  FixedBoundary.i class definition.
 */
//---------------------------------------------------------------------------//

#ifndef detran_FIXEDBOUNDARY_I_HH_
#define detran_FIXEDBOUNDARY_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
template <class D>
FixedBoundary<D>::FixedBoundary(SP_boundary boundary,
                                const size_t side,
                                SP_input input,
                                SP_mesh mesh,
                                SP_quadrature quadrature)
  : Base(boundary, side, input, mesh, quadrature)
  , d_number_groups(d_input->template get<int>("number_groups"))
  , d_psi(d_number_groups,
          vec2_bflux(d_quadrature->number_octants()/2,
                     vec1_bflux(d_quadrature->number_angles_octant(),
                                (*d_boundary)(side, 0, 0, 0))))
{
  /* ... */
}


//---------------------------------------------------------------------------//
template <class D>
inline void FixedBoundary<D>::set(const size_t g)
{
  for (size_t io = 0; io < d_quadrature->incident_octant(d_side).size(); ++io)
  {
    size_t o = d_quadrature->incident_octant(d_side)[io];
    for (size_t a = 0; a < d_quadrature->number_angles_octant(); ++a)
    {
      (*d_boundary)(d_side, o, a, g) = d_psi[g][io][a];
    }
  }
}

//---------------------------------------------------------------------------//
template <class D>
inline const typename FixedBoundary<D>::bf_type&
FixedBoundary<D>::operator()(const size_t o,
                             const size_t a,
                             const size_t g) const
{
//  std::cout << d_psi.size() << " "
//            << d_psi[0].size() << std::endl;
  Require(g < d_number_groups);
  Require(o < d_psi[g].size());
  Require(a < d_psi[g][o].size());
  return d_psi[g][o][a];
}

//---------------------------------------------------------------------------//
template <class D>
inline typename FixedBoundary<D>::bf_type&
FixedBoundary<D>::operator()(const size_t o,
                             const size_t a,
                             const size_t g)
{
  // Cast away return type
  return const_cast<typename FixedBoundary<D>::bf_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const FixedBoundary<D>&>(*this)(o, a, g)
  );
}


//---------------------------------------------------------------------------//
template <class D>
void FixedBoundary<D>::clear()
{
  size_t dim = d_side / 2;
  size_t dims[3][2] = { {1, 2}, {0, 2}, {0, 1} };
  for (size_t g = 0; g < d_psi.size(); ++g)
  {
    for (size_t o = 0; o < d_psi[g].size(); ++o)
    {
      for (size_t a = 0; a < d_psi[g][o].size(); ++a)
      {
        for (size_t i = 0; i < d_mesh->number_cells(dims[dim][0]); ++i)
        {
          for (size_t j = 0; j < d_mesh->number_cells(dims[dim][1]); ++j)
          {
            BoundaryValue<D>::value(d_psi[g][o][a], i, j) = 0.0;
          }
        }
      }
    }
  }
}

} // end namespace detran

#endif /* detran_FIXEDBOUNDARY_I_HH_ */
