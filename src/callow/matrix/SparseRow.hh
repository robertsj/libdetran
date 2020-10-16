//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  SparseRow.hh
 *  @brief Sparse row class definition.
 *  @note  Copyright(C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#ifndef callow_SPARSEROW_HH_
#define callow_SPARSEROW_HH_

#include <ostream>
#include <vector>
#include <string>
#include <utility>
#include "Matrix.hh"
#include <iostream>

namespace callow
{

/**
 *  @class SparseRow
 *  @brief Single, sparse row vector.
 *
 *  This is a helper class for ILUT, but it's becoming a tutorial
 *  for some C++ 11+ features.
 */
class CALLOW_EXPORT SparseRow
{

public:

  enum insert_type
  {
    INSERT, ADD, END_INSERT_TYPE
  };

  typedef std::pair<int, double> element_t;
  typedef std::vector<element_t> vec_element_t;
  typedef vec_element_t::iterator iterator;
  typedef vec_element_t::const_iterator citerator;
  typedef std::pair<iterator, iterator> range_t;
  typedef std::pair<citerator, citerator> crange_t;

  //--------------------------------------------------------------------------//
  // CONSTRUCTOR & DESTRUCTOR
  //--------------------------------------------------------------------------//

  // construction with sizing.
  SparseRow(const int n);
  // construction via copy of another row.
  SparseRow(const SparseRow &r);
  // construction via copy of CRS matrix's ith row.
  SparseRow(const Matrix &A, const int i);
  // construction via arbitrary array of length n
  SparseRow(const double *a, const int n);

  //--------------------------------------------------------------------------//
  // PUBLIC FUNCTIONS
  //--------------------------------------------------------------------------//

  /// add one value to actual index i.  return false if element deleted.
  bool insert(int i, double v, const int type = INSERT);

  /// this[a:b] = this[a:b]*type + v*r[a:b].  return false if any deletions.
  bool insert(int a, int b, const SparseRow &r, double v = 1.0, const int t = INSERT);

  /// delete an element using cardinal
  void delete_value(const int p);

  /// finalize with a culling of zeros and sort
  bool update(const double threshold = 0.0);

  /// get logical size
  int size() const {return d_n;}

  /// Get location of index
  int find(const int i) const;

  /// value at an index
  const double& operator[](const int i) const;
  double& operator[](const int i);

  /// get all elements
  const std::vector<element_t>& elements() const {return d_elements;}
  std::vector<element_t>& elements() {return d_elements;}

  /// number of nonzeros
  int number_nonzeros() const { return d_elements.size(); }

  /// pretty print
  void display(const std::string s = "") const;

  // sparse set to zero
  void clear();

  // Checks validity of row order
  bool verify() const;

  // Return a view with column values that are in [low, high]
 // SparseRow operator()(const int low, const int high);


  static bool compare_v(const element_t& e0, const element_t& e1)
  {
    return std::abs(e0.second) > std::abs(e1.second);
  }

  static bool compare_c(const element_t& e0, const element_t& e1)
  {
    return   e0.first < e1.first;
  }

  iterator begin() { return d_elements.begin(); }
  citerator cbegin() const { return d_elements.cbegin(); }
  iterator end() { return d_elements.end(); }
  citerator cend() const { return d_elements.cend(); }

  // return begin and end iterators for columns in [a, b]
  range_t range(const int a, const int b);
  crange_t range(const int a, const int b) const;


  /// Delete elements

private:

  //--------------------------------------------------------------------------//
  // DATA
  //--------------------------------------------------------------------------//

  /// Logical size
  int d_n;
  /// Non zero locations
  vec_element_t d_elements;

};

} // end namespace callow

#endif /* callow_SPARSEROW_HH_ */

//----------------------------------------------------------------------------//
//              end of SPARSEROW.cc
//----------------------------------------------------------------------------//
