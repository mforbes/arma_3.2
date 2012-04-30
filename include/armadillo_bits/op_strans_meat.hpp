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


//! \addtogroup op_strans
//! @{



//! for tiny square matrices (size <= 4x4)
template<typename eT>
inline
void
op_strans::apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A)
  {
  const eT* Am   = A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



//! Immediate transpose of a dense matrix
template<typename eT>
inline
void
op_strans::apply_noalias(Mat<eT>& out, const Mat<eT>& A)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_cols = A.n_cols;
  const uword A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    arrayops::copy( out.memptr(), A.mem, A.n_elem );
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans::apply_noalias_tinysq(out, A);
      }
    else
      {
      for(uword k=0; k < A_n_cols; ++k)
        {
        uword i, j;
        
        const eT* colptr = A.colptr(k);
        
        for(i=0, j=1; j < A_n_rows; i+=2, j+=2)
          {
          const eT tmp_i = colptr[i];
          const eT tmp_j = colptr[j];
          
          out.at(k, i) = tmp_i;
          out.at(k, j) = tmp_j;
          }
        
        if(i < A_n_rows)
          {
          out.at(k, i) = colptr[i];
          }
        }
      }
    }
  }



template<typename eT>
inline
void
op_strans::apply(Mat<eT>& out, const Mat<eT>& A, const bool is_alias)
  {
  arma_extra_debug_sigprint();
  
  if(is_alias == false)
    {
    op_strans::apply_noalias(out, A);
    }
  else
    {
    const uword n_rows = out.n_rows;
    const uword n_cols = out.n_cols;
    
    if(n_rows == n_cols)
      {
      arma_extra_debug_print("op_strans::apply(): doing in-place transpose of a square matrix");
      
      const uword N = n_rows;
      
      for(uword k=0; k < N; ++k)
        {
        eT* colptr = out.colptr(k);
        
        uword i,j;
        
        for(i=(k+1), j=(k+2); j < N; i+=2, j+=2)
          {
          std::swap(out.at(k,i), colptr[i]);
          std::swap(out.at(k,j), colptr[j]);
          }
        
        if(i < N)
          {
          std::swap(out.at(k,i), colptr[i]);
          }
        }
      }
    else
      {
      Mat<eT> tmp;
      op_strans::apply_noalias(tmp, A);
      
      out.steal_mem(tmp);
      }
    }
  }



template<typename T1>
inline
void
op_strans::apply(Mat<typename T1::elem_type>& out, const Op<T1,op_strans>& in)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const unwrap<T1>   tmp(in.m);
  const Mat<eT>& A = tmp.M;
  
  op_strans::apply(out, A, tmp.is_alias(out));
  }



//
// op_strans2



//! for tiny square matrices (size <= 4x4)
template<typename eT>
inline
void
op_strans2::apply_noalias_tinysq(Mat<eT>& out, const Mat<eT>& A, const eT val)
  {
  const eT* Am   = A.memptr();
        eT* outm = out.memptr();
  
  switch(A.n_rows)
    {
    case 1:
      {
      outm[0] = val * Am[0];
      }
      break;
      
    case 2:
      {
      outm[pos<false,0,0>::n2] = val * Am[pos<true,0,0>::n2];
      outm[pos<false,1,0>::n2] = val * Am[pos<true,1,0>::n2];
      
      outm[pos<false,0,1>::n2] = val * Am[pos<true,0,1>::n2];
      outm[pos<false,1,1>::n2] = val * Am[pos<true,1,1>::n2];
      }
      break;
    
    case 3:
      {
      outm[pos<false,0,0>::n3] = val * Am[pos<true,0,0>::n3];
      outm[pos<false,1,0>::n3] = val * Am[pos<true,1,0>::n3];
      outm[pos<false,2,0>::n3] = val * Am[pos<true,2,0>::n3];
      
      outm[pos<false,0,1>::n3] = val * Am[pos<true,0,1>::n3];
      outm[pos<false,1,1>::n3] = val * Am[pos<true,1,1>::n3];
      outm[pos<false,2,1>::n3] = val * Am[pos<true,2,1>::n3];
      
      outm[pos<false,0,2>::n3] = val * Am[pos<true,0,2>::n3];
      outm[pos<false,1,2>::n3] = val * Am[pos<true,1,2>::n3];
      outm[pos<false,2,2>::n3] = val * Am[pos<true,2,2>::n3];
      }
      break;
    
    case 4:
      {
      outm[pos<false,0,0>::n4] = val * Am[pos<true,0,0>::n4];
      outm[pos<false,1,0>::n4] = val * Am[pos<true,1,0>::n4];
      outm[pos<false,2,0>::n4] = val * Am[pos<true,2,0>::n4];
      outm[pos<false,3,0>::n4] = val * Am[pos<true,3,0>::n4];
      
      outm[pos<false,0,1>::n4] = val * Am[pos<true,0,1>::n4];
      outm[pos<false,1,1>::n4] = val * Am[pos<true,1,1>::n4];
      outm[pos<false,2,1>::n4] = val * Am[pos<true,2,1>::n4];
      outm[pos<false,3,1>::n4] = val * Am[pos<true,3,1>::n4];
      
      outm[pos<false,0,2>::n4] = val * Am[pos<true,0,2>::n4];
      outm[pos<false,1,2>::n4] = val * Am[pos<true,1,2>::n4];
      outm[pos<false,2,2>::n4] = val * Am[pos<true,2,2>::n4];
      outm[pos<false,3,2>::n4] = val * Am[pos<true,3,2>::n4];
      
      outm[pos<false,0,3>::n4] = val * Am[pos<true,0,3>::n4];
      outm[pos<false,1,3>::n4] = val * Am[pos<true,1,3>::n4];
      outm[pos<false,2,3>::n4] = val * Am[pos<true,2,3>::n4];
      outm[pos<false,3,3>::n4] = val * Am[pos<true,3,3>::n4];
      }
      break;
    
    default:
      ;
    }
  
  }



//! Immediate transpose of a dense matrix
template<typename eT>
inline
void
op_strans2::apply_noalias(Mat<eT>& out, const Mat<eT>& A, const eT val)
  {
  arma_extra_debug_sigprint();
  
  const uword A_n_cols = A.n_cols;
  const uword A_n_rows = A.n_rows;
  
  out.set_size(A_n_cols, A_n_rows);
  
  if( (A_n_cols == 1) || (A_n_rows == 1) )
    {
    const uword N = A.n_elem;
    
    const eT* A_mem   = A.memptr();
          eT* out_mem = out.memptr();
    
    uword i,j;
    for(i=0, j=1; j < N; i+=2, j+=2)
      {
      const eT tmp_i = A_mem[i];
      const eT tmp_j = A_mem[j];
      
      out_mem[i] = val * tmp_i;
      out_mem[j] = val * tmp_j;
      }
    
    if(i < N)
      {
      out_mem[i] = val * A_mem[i];
      }
    }
  else
    {
    if( (A_n_rows <= 4) && (A_n_rows == A_n_cols) )
      {
      op_strans2::apply_noalias_tinysq(out, A, val);
      }
    else
      {
      for(uword k=0; k < A_n_cols; ++k)
        {
        uword i, j;
        
        const eT* colptr = A.colptr(k);
        
        for(i=0, j=1; j < A_n_rows; i+=2, j+=2)
          {
          const eT tmp_i = colptr[i];
          const eT tmp_j = colptr[j];
          
          out.at(k, i) = val * tmp_i;
          out.at(k, j) = val * tmp_j;
          }
        
        if(i < A_n_rows)
          {
          out.at(k, i) = val * colptr[i];
          }
        }
      }
    }
  }



template<typename eT>
inline
void
op_strans2::apply(Mat<eT>& out, const Mat<eT>& A, const eT val, const bool is_alias)
  {
  arma_extra_debug_sigprint();
  
  if(is_alias == false)
    {
    op_strans2::apply_noalias(out, A, val);
    }
  else
    {
    const uword n_rows = out.n_rows;
    const uword n_cols = out.n_cols;
    
    if(n_rows == n_cols)
      {
      arma_extra_debug_print("op_strans2::apply(): doing in-place transpose of a square matrix");
      
      const uword N = n_rows;
      
      // TODO: do multiplication while swapping
      
      for(uword k=0; k < N; ++k)
        {
        eT* colptr = out.colptr(k);
        
        uword i,j;
        
        for(i=(k+1), j=(k+2); j < N; i+=2, j+=2)
          {
          std::swap(out.at(k,i), colptr[i]);
          std::swap(out.at(k,j), colptr[j]);
          }
        
        if(i < N)
          {
          std::swap(out.at(k,i), colptr[i]);
          }
        }
        
      arrayops::inplace_mul( out.memptr(), val, out.n_elem );
      }
    else
      {
      Mat<eT> tmp;
      op_strans2::apply_noalias(tmp, A, val);
      
      out.steal_mem(tmp);
      }
    }
  }



//! @}
