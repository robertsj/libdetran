//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   Definitions.hh
 * \author Jeremy Roberts
 * \brief  Useful type definitions, etc.
 */
//---------------------------------------------------------------------------//

#ifndef DETRAN_UTILS_DEFINITIONS_HH_
#define DETRAN_UTILS_DEFINITIONS_HH_

#include <vector>

namespace detran_utilities
{

typedef unsigned int size_t;
typedef unsigned int u_int;

// Templated standard vec

template <typename T>
struct vec1_T {typedef std::vector<T> value_type;};
template <typename T>
struct vec2_T {typedef std::vector<typename vec1_T<T>::value_type> value_type;};
template <typename T>
struct vec3_T {typedef std::vector<typename vec2_T<T>::value_type> value_type;};
template <typename T>
struct vec4_T {typedef std::vector<typename vec3_T<T>::value_type> value_type;};

// Instantiations

//typedef vec1_T<int>::value_type      vec_int;
//typedef vec2_T<int>::value_type     vec2_int;
//typedef vec3_T<int>::value_type     vec3_int;
//typedef vec4_T<int>::value_type     vec4_int;
//
//typedef vec1_T<double>::value_type      vec_dbl;
//typedef vec2_T<double>::value_type     vec2_dbl;
//typedef vec3_T<double>::value_type     vec3_dbl;
//typedef vec4_T<double>::value_type     vec4_dbl;

typedef std::vector<int> vec_int;
typedef std::vector<std::vector<int> > vec2_int;
typedef std::vector<std::vector<std::vector<int> > > vec3_int;
typedef std::vector<double> vec_dbl;
typedef std::vector<std::vector<double> > vec2_dbl;
typedef std::vector<std::vector<std::vector<double> > > vec3_dbl;


typedef std::vector<size_t> vec_size_t;
typedef std::vector<std::vector<size_t> > vec2_size_t;
typedef std::vector<std::vector<std::vector<size_t> > > vec3_size_t;
typedef std::vector<std::vector<std::vector<std::vector<size_t> > > > vec4_size_t;
//typedef vec1_T<size_t>::value_type      vec_size_t;
//typedef vec2_T<size_t>::value_type     vec2_size_t;
//typedef vec3_T<size_t>::value_type     vec3_size_t;
//typedef vec4_T<size_t>::value_type     vec4_size_t;



} // end namespace detran_utilities


#endif // DETRAN_UTILS_DEFINITIONS_HH_
