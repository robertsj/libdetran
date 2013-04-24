//----------------------------------*-C++-*----------------------------------//
/**
 *  @file   State.i.hh
 *  @author robertsj
 *  @date   Sep 7, 2012
 *  @brief  State inline member definitions
 */
//---------------------------------------------------------------------------//

#ifndef detran_STATE_I_HH_
#define detran_STATE_I_HH_

namespace detran
{

//---------------------------------------------------------------------------//
inline const State::moments_type&
State::phi(const size_t g) const
{
  Require(g < d_number_groups);
  return d_moments[g];
}

//---------------------------------------------------------------------------//
inline State::moments_type&
State::phi(const size_t g)
{
  // Cast away return type
  return const_cast<moments_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const State*>(this)->phi(g)
  );
}

//---------------------------------------------------------------------------//
inline const State::group_moments_type&
State::all_phi() const
{
  return d_moments;
}

//---------------------------------------------------------------------------//
inline State::group_moments_type&
State::all_phi()
{
  // Cast away return type
  return const_cast<group_moments_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const State*>(this)->all_phi()
  );

}

//---------------------------------------------------------------------------//
inline void State::set_moments(const size_t g, std::vector<double>& f)
{
  d_moments[g] = f;
}

//---------------------------------------------------------------------------//
inline const State::angular_flux_type&
State::psi(const size_t g, const size_t o, const size_t a) const
{
  Require(d_store_angular_flux);
  Require(d_angular_flux.size() > 0);
  Require(o < d_quadrature->number_octants());
  Require(a < d_quadrature->number_angles_octant());
  Require(g < d_number_groups);
  int angle = d_quadrature->index(o, a);
  return d_angular_flux[g][angle];
}

//---------------------------------------------------------------------------//
inline State::angular_flux_type&
State::psi(const size_t g, const size_t o, const size_t a)
{
  // Cast away return type
  return const_cast<angular_flux_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const State*>(this)->psi(g, o, a)
  );
}

//---------------------------------------------------------------------------//
inline const State::moments_type&
State::current(const size_t g) const
{
  Require(d_store_current);
  Require(g < d_number_groups);
  return d_current[g];
}

//---------------------------------------------------------------------------//
inline State::moments_type&
State::current(const size_t g)
{
  // Cast away return type
  return const_cast<moments_type&>
  (
    // Add const to *this's type and call const version
    static_cast<const State*>(this)->current(g)
  );
}

} // end namespace detran


#endif /* detran_STATE_I_HH_ */

//---------------------------------------------------------------------------//
//              end of State.i.hh
//---------------------------------------------------------------------------//
