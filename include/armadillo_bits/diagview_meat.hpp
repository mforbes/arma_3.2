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


//! \addtogroup diagview
//! @{


template<typename eT>
inline
diagview<eT>::~diagview()
  {
  arma_extra_debug_sigprint();
  }


template<typename eT>
arma_inline
diagview<eT>::diagview(const Mat<eT>& in_m, const uword in_row_offset, const uword in_col_offset, const uword in_len)
  : m(in_m)
  , row_offset(in_row_offset)
  , col_offset(in_col_offset)
  , n_rows(in_len)
  , n_elem(in_len)
  {
  arma_extra_debug_sigprint();
  }



//! set a diagonal of our matrix using a diagonal from a foreign matrix
template<typename eT>
inline
void
diagview<eT>::operator= (const diagview<eT>& x)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  arma_debug_check( (d.n_elem != x.n_elem), "diagview: diagonals have incompatible lengths");
  
        Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  const Mat<eT>& x_m = x.m;
  
  if(&d_m != &x_m)
    {
    const uword d_n_elem     = d.n_elem;
    const uword d_row_offset = d.row_offset;
    const uword d_col_offset = d.col_offset;
    
    const uword x_row_offset = x.row_offset;
    const uword x_col_offset = x.col_offset;
    
    uword i,j;
    for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
      {
      const eT tmp_i = x_m.at(i + x_row_offset, i + x_col_offset);
      const eT tmp_j = x_m.at(j + x_row_offset, j + x_col_offset);
      
      d_m.at(i + d_row_offset, i + d_col_offset) = tmp_i;
      d_m.at(j + d_row_offset, j + d_col_offset) = tmp_j;
      }
    
    if(i < d_n_elem)
      {
      d_m.at(i + d_row_offset, i + d_col_offset) = x_m.at(i + x_row_offset, i + x_col_offset);
      }
    }
  else
    {
    const Mat<eT> tmp = x;
    
    (*this).operator=(tmp);
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator+=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = const_cast< Mat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) += val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator-=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = const_cast< Mat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) -= val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator*=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = const_cast< Mat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) *= val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::operator/=(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& t_m = const_cast< Mat<eT>& >(m);
  
  const uword t_n_elem     = n_elem;
  const uword t_row_offset = row_offset;
  const uword t_col_offset = col_offset;
  
  for(uword i=0; i<t_n_elem; ++i)
    {
    t_m.at( i + t_row_offset,  i + t_col_offset) /= val;
    }
  }



//! set a diagonal of our matrix using data from a foreign object
template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator= (const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  const unwrap_check<T1> tmp( o.get_ref(), d.m );
  const Mat<eT>& x = tmp.M;
  
  arma_debug_check
    (
    ( (d.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    d_m.at( i + d_row_offset,  i + d_col_offset) = tmp_i;
    d_m.at( j + d_row_offset,  j + d_col_offset) = tmp_j;
    }
  
  if(i < d_n_elem)
    {
    d_m.at( i + d_row_offset,  i + d_col_offset) = x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator+=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  const unwrap_check<T1> tmp( o.get_ref(), d.m );
  const Mat<eT>& x = tmp.M;
  
  arma_debug_check
    (
    ( (d.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    d_m.at( i + d_row_offset,  i + d_col_offset) += tmp_i;
    d_m.at( j + d_row_offset,  j + d_col_offset) += tmp_j;
    }
  
  if(i < d_n_elem)
    {
    d_m.at( i + d_row_offset,  i + d_col_offset) += x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator-=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  const unwrap_check<T1> tmp( o.get_ref(), d.m );
  const Mat<eT>& x = tmp.M;
  
  arma_debug_check
    (
    ( (d.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    d_m.at( i + d_row_offset,  i + d_col_offset) -= tmp_i;
    d_m.at( j + d_row_offset,  j + d_col_offset) -= tmp_j;
    }
  
  if(i < d_n_elem)
    {
    d_m.at( i + d_row_offset,  i + d_col_offset) -= x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator%=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  const unwrap_check<T1> tmp( o.get_ref(), d.m );
  const Mat<eT>& x = tmp.M;
  
  arma_debug_check
    (
    ( (d.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    d_m.at( i + d_row_offset,  i + d_col_offset) *= tmp_i;
    d_m.at( j + d_row_offset,  j + d_col_offset) *= tmp_j;
    }
  
  if(i < d_n_elem)
    {
    d_m.at( i + d_row_offset,  i + d_col_offset) *= x_mem[i];
    }
  }



template<typename eT>
template<typename T1>
inline
void
diagview<eT>::operator/=(const Base<eT,T1>& o)
  {
  arma_extra_debug_sigprint();
  
  diagview<eT>& d = *this;
  
  const unwrap_check<T1> tmp( o.get_ref(), d.m );
  const Mat<eT>& x = tmp.M;
  
  arma_debug_check
    (
    ( (d.n_elem != x.n_elem) || (x.is_vec() == false) ),
    "diagview: given object has incompatible size"
    );
  
  Mat<eT>& d_m = const_cast< Mat<eT>& >(d.m);
  
  const uword d_n_elem     = d.n_elem;
  const uword d_row_offset = d.row_offset;
  const uword d_col_offset = d.col_offset;
  
  const eT* x_mem = x.memptr();
  
  uword i,j;
  for(i=0, j=1; j < d_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = x_mem[i];
    const eT tmp_j = x_mem[j];
    
    d_m.at( i + d_row_offset,  i + d_col_offset) /= tmp_i;
    d_m.at( j + d_row_offset,  j + d_col_offset) /= tmp_j;
    }
  
  if(i < d_n_elem)
    {
    d_m.at( i + d_row_offset,  i + d_col_offset) /= x_mem[i];
    }
  }



//! extract a diagonal and store it as a column vector
template<typename eT>
inline
void
diagview<eT>::extract(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  // NOTE: we're assuming that the matrix has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Mat contructor or operator=()
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] = tmp_i;
    out_mem[j] = tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] = in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X += Y.diag()
template<typename eT>
inline
void
diagview<eT>::plus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "addition");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] += tmp_i;
    out_mem[j] += tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] += in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X -= Y.diag()
template<typename eT>
inline
void
diagview<eT>::minus_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "subtraction");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] -= tmp_i;
    out_mem[j] -= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] -= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X %= Y.diag()
template<typename eT>
inline
void
diagview<eT>::schur_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "element-wise multiplication");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] *= tmp_i;
    out_mem[j] *= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] *= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



//! X /= Y.diag()
template<typename eT>
inline
void
diagview<eT>::div_inplace(Mat<eT>& out, const diagview<eT>& in)
  {
  arma_extra_debug_sigprint();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, in.n_rows, in.n_cols, "element-wise division");
  
  const Mat<eT>& in_m = in.m;
  
  const uword in_n_elem     = in.n_elem;
  const uword in_row_offset = in.row_offset;
  const uword in_col_offset = in.col_offset;
  
  eT* out_mem = out.memptr();
  
  uword i,j;
  for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
    {
    const eT tmp_i = in_m.at( i + in_row_offset, i + in_col_offset );
    const eT tmp_j = in_m.at( j + in_row_offset, j + in_col_offset );
    
    out_mem[i] /= tmp_i;
    out_mem[j] /= tmp_j;
    }
  
  if(i < in_n_elem)
    {
    out_mem[i] /= in_m.at( i + in_row_offset, i + in_col_offset );
    }
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator[](const uword ii)
  {
  return (const_cast< Mat<eT>& >(m)).at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator[](const uword ii) const
  {
  return m.at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::at(const uword ii)
  {
  return (const_cast< Mat<eT>& >(m)).at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::at(const uword ii) const
  {
  return m.at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const uword ii)
  {
  arma_debug_check( (ii >= n_elem), "diagview::operator(): out of bounds" );
  
  return (const_cast< Mat<eT>& >(m)).at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const uword ii) const
  {
  arma_debug_check( (ii >= n_elem), "diagview::operator(): out of bounds" );
  
  return m.at(ii+row_offset, ii+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::at(const uword row, const uword)
  {
  return (const_cast< Mat<eT>& >(m)).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::at(const uword row, const uword) const
  {
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT&
diagview<eT>::operator()(const uword row, const uword col)
  {
  arma_debug_check( ((row >= n_elem) || (col > 0)), "diagview::operator(): out of bounds" );
  
  return (const_cast< Mat<eT>& >(m)).at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
eT
diagview<eT>::operator()(const uword row, const uword col) const
  {
  arma_debug_check( ((row >= n_elem) || (col > 0)), "diagview::operator(): out of bounds" );
  
  return m.at(row+row_offset, row+col_offset);
  }



template<typename eT>
arma_inline
const Op<diagview<eT>,op_htrans>
diagview<eT>::t() const
  {
  return Op<diagview<eT>,op_htrans>(*this);
  }



template<typename eT>
arma_inline
const Op<diagview<eT>,op_htrans>
diagview<eT>::ht() const
  {
  return Op<diagview<eT>,op_htrans>(*this);
  }



template<typename eT>
arma_inline
const Op<diagview<eT>,op_strans>
diagview<eT>::st() const
  {
  return Op<diagview<eT>,op_strans>(*this);
  }



template<typename eT>
inline
void
diagview<eT>::fill(const eT val)
  {
  arma_extra_debug_sigprint();
  
  Mat<eT>& x = const_cast< Mat<eT>& >(m);
  
  for(uword ii=0; ii < n_elem; ++ii)
    {
    x.at(ii+row_offset, ii+col_offset) = val;
    }
  }



template<typename eT>
inline
void
diagview<eT>::zeros()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(0));
  }



template<typename eT>
inline
void
diagview<eT>::ones()
  {
  arma_extra_debug_sigprint();
  
  (*this).fill(eT(1));
  }



//! @}
