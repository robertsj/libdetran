//----------------------------------*-C++-*----------------------------------//
/**
 *  @file  DimensionTraits.hh
 *  @brief DimensionTraits class definition
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//---------------------------------------------------------------------------//

#ifndef detran_DIMENSIONTRAITS_HH_
#define detran_DIMENSIONTRAITS_HH_

namespace detran
{

//---------------------------------------------------------------------------//
/**
 *  @class DimensionTraits
 *  @brief Allows templates using dimension without <int N> format.
 *
 *  There may eventually be further uses for the class.
 */
//---------------------------------------------------------------------------//
template <int N>
class DimensionTraits
{
public:
    typedef int value_type;
    static const value_type dimension = N;
};

// instantiations
typedef DimensionTraits<1> _1D;
typedef DimensionTraits<1> _1DCYL; // a reminder to explore other options
typedef DimensionTraits<1> _1DSPH;
typedef DimensionTraits<2> _2D;
typedef DimensionTraits<2> _2DHEX; // ...
typedef DimensionTraits<3> _3D;

} // end namespace detran

#endif /* detran_DIMENSIONTRAITS_HH_ */

//----------------------------------------------------------------------------//
//              end of file DimensionTraits.hh
//----------------------------------------------------------------------------//
