// Copyright (C) 2012 NICTA (www.nicta.com.au)
// Copyright (C) 2012 Conrad Sanderson
// Copyright (C) 2012 Arnold Wiliem
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)



//! \addtogroup op_unique
//! @{



// TODO: add an efficient implementation for complex numbers

template<typename T1>
inline
void
op_unique::apply(Mat<typename T1::elem_type>& out, const Op<T1, op_unique>& X)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const Proxy<T1> P(X.m);
  
  const uword in_n_rows = P.get_n_rows();
  const uword in_n_cols = P.get_n_cols();
  const uword in_n_elem = P.get_n_elem();
  
  
  if(in_n_elem <= 1)
    {
    if(in_n_elem == 1)
      {
      const eT tmp = P[0];
      
      out.set_size(in_n_rows, in_n_cols);
      
      out[0] = tmp;
      }
    else
      {
      out.set_size(in_n_rows, in_n_cols);
      }
    
    return;
    }
  
  
  std::vector<eT> vec(in_n_elem);
  
  
  if(Proxy<T1>::prefer_at_accessor == false)
    {
    typename Proxy<T1>::ea_type Pea = P.get_ea();
    
    uword i,j;
    for(i=0, j=1; j < in_n_elem; i+=2, j+=2)
      {
      const eT tmp_i = Pea[i];
      const eT tmp_j = Pea[j];
      
      vec[i] = tmp_i;
      vec[j] = tmp_j;
      }
    
    if(i < in_n_elem)
      {
      vec[i] = Pea[i];
      }
    }
  else
    {
    uword i = 0;
    
    for(uword col=0; col < in_n_cols; ++col)
    for(uword row=0; row < in_n_rows; ++row, ++i)
      {
      vec[i] = P.at(row,col);
      }
    }
  
  std::sort( vec.begin(), vec.end() );
  
  uword N_unique = 1;
  
  for(uword i=1; i < in_n_elem; ++i)
    {
    const eT a = vec[i-1];
    const eT b = vec[i  ];
    
    const eT diff = a - b;
    
    if(diff != eT(0)) { ++N_unique; }
    }
  
  uword out_n_rows;
  uword out_n_cols;
  
  if( (in_n_rows == 1) || (in_n_cols == 1) )
    {
    if(in_n_rows == 1)
      {
      out_n_rows = 1;
      out_n_cols = N_unique;
      }
    else
      {
      out_n_rows = N_unique;
      out_n_cols = 1;
      }
    }
  else
    {
    out_n_rows = N_unique;
    out_n_cols = 1;
    }
  
  // we don't need to worry about aliasing at this stage, as all the data is stored in vec
  out.set_size(out_n_rows, out_n_cols);
  
  eT* out_mem = out.memptr();
  
  if(in_n_elem > 0) { out_mem[0] = vec[0]; }
  
  N_unique = 1;
  
  for(uword i=1; i < in_n_elem; ++i)
    {
    const eT a = vec[i-1];
    const eT b = vec[i  ];
    
    const eT diff = a - b;
    
    if(diff != eT(0))
      {
      out_mem[N_unique] = b;
      ++N_unique;
      }
    }
  
  }



//! @}
