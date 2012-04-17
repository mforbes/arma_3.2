// Copyright (C) 2008-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup op_diagmat
//! @{



template<typename T1>
inline
void
op_diagmat::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_diagmat>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.m);
  
  const uword n_rows = P.get_n_rows();
  const uword n_cols = P.get_n_cols();
  
  if(P.is_alias(out) == false)
    {
    if( (n_rows == 1) || (n_cols == 1) )    // generate a diagonal matrix out of a vector
      {
      if(n_rows == 1)
        {
        out.zeros(n_cols, n_cols);
        
        for(uword i=0; i < n_cols; ++i)
          {
          out.at(i,i) = (Proxy<T1>::prefer_at_accessor == false) ? P[i] : P.at(0,i);
          }
        }
      else
      if(n_cols == 1)
        {
        out.zeros(n_rows, n_rows);
        
        for(uword i=0; i < n_rows; ++i)
          {
          out.at(i,i) = (Proxy<T1>::prefer_at_accessor == false) ? P[i] : P.at(i,0);
          }
        }
      }
    else   // generate a diagonal matrix out of a matrix
      {
      arma_debug_check( (n_rows != n_cols), "diagmat(): given matrix is not square" );
      
      out.zeros(n_rows, n_rows);
      
      for(uword i=0; i < n_rows; ++i)
        {
        out.at(i,i) = P.at(i,i);
        }
      }
    }
  else   // we have aliasing
    {
    if( (n_rows == 1) || (n_cols == 1) )   // generate a diagonal matrix out of a vector
      {
      podarray<eT> tmp;
      
      eT* tmp_mem;
      
      if(n_rows == 1)
        {
        tmp.set_size(n_cols);
        
        tmp_mem = tmp.memptr();
        
        for(uword i=0; i < n_cols; ++i)
          {
          tmp_mem[i] = (Proxy<T1>::prefer_at_accessor == false) ? P[i] : P.at(0,i);
          }
        }
      else
      if(n_cols == 1)
        {
        tmp.set_size(n_rows);
        
        tmp_mem = tmp.memptr();
        
        for(uword i=0; i < n_rows; ++i)
          {
          tmp_mem[i] = (Proxy<T1>::prefer_at_accessor == false) ? P[i] : P.at(i,0);
          }
        }
      
      
      const uword n_elem = tmp.n_elem;
      
      out.zeros(n_elem, n_elem);
      
      for(uword i=0; i < n_elem; ++i)
        {
        out.at(i,i) = tmp_mem[i];
        }
      }
    else   // generate a diagonal matrix out of a matrix
      {
      arma_debug_check( (n_rows != n_cols), "diagmat(): given matrix is not square" );
      
      for(uword i=0; i < n_rows; ++i)
        {
        eT* colptr = out.colptr(i);
        
        // clear above the diagonal
        arrayops::inplace_set(colptr, eT(0), i);
        
        // clear below the diagonal
        arrayops::inplace_set(colptr+(i+1), eT(0), n_rows-1-i);
        }
      }
    }
  }



//! @}
