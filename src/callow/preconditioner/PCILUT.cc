//----------------------------------*-C++-*-----------------------------------//
/**
 *  @file  PCILUT.cc
 *  @brief PCILUT member definitions
 *  @note  Copyright (C) 2012-2013 Jeremy Roberts
 */
//----------------------------------------------------------------------------//

#include "PCILUT.hh"
#include "SparseRow.hh"

#include <cmath>
#include <algorithm>
#include <iostream>
using std::cout;
using std::endl;

namespace callow
{


  bool drop_t(SparseRow  &r, int idx, double t)
  {
    if (std::abs(r[idx]) < 0.0)
    {
      r.delete_value(idx);
      return false;
    }
    return true;
  }

  bool drop_p(SparseRow  &r, int d, int p, double t)
  {
    // first drop small values
    for (int idx = 0; idx < r.number_nonzeros(); ++idx)
    {
      drop_t(r, idx, t);
    }
    // then keep the p largest on either side of d
    int cL = 0;
    int idx_d = 0;
    for (int idx = 0; idx < r.number_nonzeros(); ++idx)
    {
      int j = r.elements()[idx].first;
      if (j == d)
        idx_d = idx;
      else if (j < d)
        cL++;
    }
    int cU = r.number_nonzeros() - 1 - cL;
    Ensure(cU >= 0);

    // copy r and sort the first and second portions
    auto q = r;
    //r.display("r");
    auto rq = q.range(0, d);
    std::sort(rq.first, rq.second, SparseRow::compare_v);
    rq = q.range(d+1, r.size());
    std::sort(rq.first, rq.second, SparseRow::compare_v);
    // iterators are invalid now.  however, we know we need to keep
    // the largest p (or drop the smallest)
    int drop = cL - p;
    if (drop > 0)
    {
    	q.elements().erase(q.begin(), q.begin()+drop);
    	cL = p;
    }
    drop = cU - p;
    if (drop>0)
    {
    	q.elements().erase(q.begin()+cL+1, q.begin()+cL+1+drop);
    	cU = p;
    }
    // sort q by column and copy to r
    std::sort(q.begin(), q.end(), SparseRow::compare_c);
    //q.display("q");
    r = q;
    return true;
  }

//----------------------------------------------------------------------------//
PCILUT::PCILUT(SP_matrix A, const size_t p, const double t)
  : Base("PCILUT")
  , d_p(p)
  , d_t(t)
{
  Require(A);
  Require(A->number_rows() == A->number_columns());
  Insist(dynamic_cast<Matrix*>(A.bp()),
    "Need an explicit matrix for use with PCILUT");
  if (d_p == 0)
    d_p = A->number_rows();

  SP_matrixfull B(A);

  d_size = A->number_columns();

  // get the csr structure
  int n = B->number_rows();

  std::vector<SparseRow> rows;
  // P = L + U

  bool keep = false;
  for (int i = 0; i < n; ++i)
  {
    SparseRow row(*B, i);
    for (int k = 0; k < i; ++k)
    {
      int p = row.find(k);
      if (p < 0)
        continue;
      SparseRow &row_U = rows[k];
      double diag_val = row_U.elements()[row_U.find(k)].second;
      Ensure(diag_val != 0);
      row.display("row befor scale");
      row.elements()[p].second /= diag_val;
      //keep = drop_t(row, p, d_t);
      //if (!keep)
      //  continue;
      //
      cout << "*************************************" << endl;
      row.display("row");
      row_U.display("P[k, :]");
      row.insert(k+1, n, row_U, -row[k], SparseRow::ADD);
      row.display("row");
      cout << " i = " << i << "  k = " << k  << " diag = " << diag_val << endl;
      cout << "*************************************" << endl;

    } // end k
    //keep = drop_p(row, i, d_p, d_t);
    //row.display("final");
    rows.push_back(row);
  }

  // copy A
  d_P = new Matrix(n, n, 2*d_p+1);
  for (int i = 0; i < n; ++i)
  {
    for (auto it = rows[i].begin(); it < rows[i].end(); ++it)
    {
      d_P->insert(i, it->first, it->second);
	}
  }
  d_P->assemble();


  cout << "...FINISHED..." << endl;
}

//----------------------------------------------------------------------------//
PCILUT::SP_preconditioner
PCILUT::Create(SP_matrix A, const size_t p, const double t)
{
  SP_preconditioner pc(new PCILUT(A, p, t));
  return pc;
}

//----------------------------------------------------------------------------//
void PCILUT::display(const std::string &name)
{
  cout << " PCILU " << name << endl;
  d_P->display();
}

//----------------------------------------------------------------------------//
void PCILUT::apply(Vector &b, Vector &x)
{
  //
}


} // end namespace callow

//----------------------------------------------------------------------------//
//              end of file PCILU0.cc
//----------------------------------------------------------------------------//
