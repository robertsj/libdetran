/*
 * Traits.hh
 *
 *  Created on: Apr 4, 2012
 *      Author: robertsj
 */

#ifndef TRAITS_HH_
#define TRAITS_HH_

namespace detran
{

//===========================================================================//
/*!
 * \class Dimension
 * \brief Allows templates using dimension without <int N> format.
 *
 * There may eventually be further uses for the class.
 *
 * Idea taken from S. Johnson's code PyTRT.
 *
 */
//===========================================================================//
template <int N>
class Dimension
{
public:
    typedef int value_type;
    static const value_type dimension = N;
};
// instantiations
typedef Dimension<1> _1D;
typedef Dimension<2> _2D;
typedef Dimension<3> _3D;

} // end namespace detran

#endif /* TRAITS_HH_ */
