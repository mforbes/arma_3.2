// Copyright (C) 2010-2012 NICTA (www.nicta.com.au)
// Copyright (C) 2010-2012 Conrad Sanderson
// 
// This file is part of the Armadillo C++ library.
// It is provided without any warranty of fitness
// for any purpose. You can redistribute this file
// and/or modify it under the terms of the GNU
// Lesser General Public License (LGPL) as published
// by the Free Software Foundation, either version 3
// of the License or (at your option) any later version.
// (see http://www.opensource.org/licenses for more info)


//! \addtogroup eglue_core
//! @{



class eglue_plus : public eglue_core<eglue_plus>
  {
  public:
  
  inline static const char* text() { return "addition"; }
  };



class eglue_minus : public eglue_core<eglue_minus>
  {
  public:
  
  inline static const char* text() { return "subtraction"; }
  };



class eglue_div : public eglue_core<eglue_div>
  {
  public:
  
  inline static const char* text() { return "element-wise division"; }
  };



class eglue_schur : public eglue_core<eglue_schur>
  {
  public:
  
  inline static const char* text() { return "element-wise multiplication"; }
  };



#undef arma_applier_1
#undef arma_applier_2
#undef arma_applier_3
#undef operatorA
#undef operatorB

#define arma_applier_1(operatorA, operatorB) \
  {\
  if( (is_fixed == true) && (n_elem <= 16) ) \
    {\
    for(uword i=0; i<n_elem; ++i)\
      {\
      out_mem[i] operatorA P1[i] operatorB P2[i];\
      }\
    }\
  else \
    {\
    uword i,j;\
    \
    for(i=0, j=1; j<n_elem; i+=2, j+=2)\
      {\
      eT tmp_i = P1[i];\
      eT tmp_j = P1[j];\
      \
      tmp_i operatorB##= P2[i];\
      tmp_j operatorB##= P2[j];\
      \
      out_mem[i] operatorA tmp_i;\
      out_mem[j] operatorA tmp_j;\
      }\
    \
    if(i < n_elem)\
      {\
      out_mem[i] operatorA P1[i] operatorB P2[i];\
      }\
    }\
  }
  


#define arma_applier_2(operatorA, operatorB) \
  {\
  if(n_rows != 1)\
    {\
    uword count = 0;\
    \
    for(uword col=0; col<n_cols; ++col)\
      {\
      uword i,j;\
      \
      for(i=0, j=1; j<n_rows; i+=2, j+=2, count+=2)\
        {\
        eT tmp_i = P1.at(i,col);\
        eT tmp_j = P1.at(j,col);\
        \
        tmp_i operatorB##= P2.at(i,col);\
        tmp_j operatorB##= P2.at(j,col);\
        \
        out_mem[count  ] operatorA tmp_i;\
        out_mem[count+1] operatorA tmp_j;\
        }\
      \
      if(i < n_rows)\
        {\
        out_mem[count] operatorA P1.at(i,col) operatorB P2.at(i,col);\
        ++count;\
        }\
      }\
    }\
  else\
    {\
    uword i,j;\
    for(i=0, j=1; j < n_cols; i+=2, j+=2)\
      {\
      eT tmp_i = P1.at(0,i);\
      eT tmp_j = P1.at(0,j);\
      \
      tmp_i operatorB##= P2.at(0,i);\
      tmp_j operatorB##= P2.at(0,j);\
      \
      out_mem[i] operatorA tmp_i;\
      out_mem[j] operatorA tmp_j;\
      }\
    \
    if(i < n_cols)\
      {\
      out_mem[i] operatorA P1.at(0,i) operatorB P2.at(0,i);\
      }\
    }\
  }



#define arma_applier_3(operatorA, operatorB) \
  {\
  uword count = 0;\
  \
  for(uword slice=0; slice<n_slices; ++slice)\
    {\
    for(uword col=0; col<n_cols; ++col)\
      {\
      uword i,j;\
      \
      for(i=0, j=1; j<n_rows; i+=2, j+=2, count+=2)\
        {\
        eT tmp_i = P1.at(i,col,slice);\
        eT tmp_j = P1.at(j,col,slice);\
        \
        tmp_i operatorB##= P2.at(i,col,slice);\
        tmp_j operatorB##= P2.at(j,col,slice);\
        \
        out_mem[count  ] operatorA tmp_i;\
        out_mem[count+1] operatorA tmp_j;\
        }\
      \
      if(i < n_rows)\
        {\
        out_mem[count] operatorA P1.at(i,col,slice) operatorB P2.at(i,col,slice);\
        ++count;\
        }\
      }\
    }\
  }



//
// matrices



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  // NOTE: we're assuming that the matrix has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Mat contructor or operator=()
  
  eT* out_mem = out.memptr();
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = (Proxy<T1>::is_fixed) || (Proxy<T2>::is_fixed);
    const uword n_elem   = x.get_n_elem();
    
    typename Proxy<T1>::ea_type P1 = x.P1.get_ea();
    typename Proxy<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(=, *); }
    }
  else
    {
    const uword n_rows = x.get_n_rows();
    const uword n_cols = x.get_n_cols();
    
    const Proxy<T1>& P1 = x.P1;
    const Proxy<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_2(=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_2(=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_2(=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_2(=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_plus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = x.get_n_rows();
  const uword n_cols = x.get_n_cols();
    
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "addition");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = (Proxy<T1>::is_fixed) || (Proxy<T2>::is_fixed);
    const uword n_elem   = x.get_n_elem();
    
    typename Proxy<T1>::ea_type P1 = x.P1.get_ea();
    typename Proxy<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(+=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(+=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(+=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(+=, *); }
    }
  else
    {
    const Proxy<T1>& P1 = x.P1;
    const Proxy<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_2(+=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_2(+=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_2(+=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_2(+=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_minus(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = x.get_n_rows();
  const uword n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "subtraction");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = (Proxy<T1>::is_fixed) || (Proxy<T2>::is_fixed);
    const uword n_elem   = x.get_n_elem();
    
    typename Proxy<T1>::ea_type P1 = x.P1.get_ea();
    typename Proxy<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(-=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(-=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(-=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(-=, *); }
    }
  else
    {
    const Proxy<T1>& P1 = x.P1;
    const Proxy<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_2(-=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_2(-=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_2(-=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_2(-=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_schur(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = x.get_n_rows();
  const uword n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise multiplication");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = (Proxy<T1>::is_fixed) || (Proxy<T2>::is_fixed);
    const uword n_elem   = x.get_n_elem();
    
    typename Proxy<T1>::ea_type P1 = x.P1.get_ea();
    typename Proxy<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(*=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(*=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(*=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(*=, *); }
    }
  else
    {
    const Proxy<T1>& P1 = x.P1;
    const Proxy<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_2(*=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_2(*=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_2(*=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_2(*=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_div(Mat<typename T1::elem_type>& out, const eGlue<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows = x.get_n_rows();
  const uword n_cols = x.get_n_cols();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, n_rows, n_cols, "element-wise division");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (Proxy<T1>::prefer_at_accessor || Proxy<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = (Proxy<T1>::is_fixed) || (Proxy<T2>::is_fixed);
    const uword n_elem   = x.get_n_elem();
    
    typename Proxy<T1>::ea_type P1 = x.P1.get_ea();
    typename Proxy<T2>::ea_type P2 = x.P2.get_ea();
  
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(/=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(/=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(/=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(/=, *); }
    }
  else
    {
    const Proxy<T1>& P1 = x.P1;
    const Proxy<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_2(*=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_2(*=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_2(*=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_2(*=, *); }
    }
  }



//
// cubes



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  typedef typename T1::elem_type eT;
  
  const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  
  // NOTE: we're assuming that the cube has already been set to the correct size and there is no aliasing;
  // size setting and alias checking is done by either the Cube contructor or operator=()
  
  
  eT* out_mem = out.memptr();
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = false;
    const uword n_elem   = out.n_elem;
    
    typename ProxyCube<T1>::ea_type P1 = x.P1.get_ea();
    typename ProxyCube<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(=, *); }
    }
  else
    {
    const uword n_rows   = x.get_n_rows();
    const uword n_cols   = x.get_n_cols();
    const uword n_slices = x.get_n_slices();
  
    const ProxyCube<T1>& P1 = x.P1;
    const ProxyCube<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_3(=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_3(=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_3(=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_3(=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_plus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows   = x.get_n_rows();
  const uword n_cols   = x.get_n_cols();
  const uword n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "addition");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = false;
    const uword n_elem   = out.n_elem;
    
    typename ProxyCube<T1>::ea_type P1 = x.P1.get_ea();
    typename ProxyCube<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(+=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(+=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(+=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(+=, *); }
    }
  else
    {
    const ProxyCube<T1>& P1 = x.P1;
    const ProxyCube<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_3(+=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_3(+=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_3(+=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_3(+=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_minus(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows   = x.get_n_rows();
  const uword n_cols   = x.get_n_cols();
  const uword n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "subtraction");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = false;
    const uword n_elem   = out.n_elem;
    
    typename ProxyCube<T1>::ea_type P1 = x.P1.get_ea();
    typename ProxyCube<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(-=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(-=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(-=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(-=, *); }
    }
  else
    {
    const ProxyCube<T1>& P1 = x.P1;
    const ProxyCube<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_3(-=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_3(-=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_3(-=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_3(-=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_schur(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows   = x.get_n_rows();
  const uword n_cols   = x.get_n_cols();
  const uword n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "element-wise multiplication");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = false;
    const uword n_elem   = out.n_elem;
    
    typename ProxyCube<T1>::ea_type P1 = x.P1.get_ea();
    typename ProxyCube<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(*=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(*=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(*=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(*=, *); }
    }
  else
    {
    const ProxyCube<T1>& P1 = x.P1;
    const ProxyCube<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_3(*=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_3(*=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_3(*=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_3(*=, *); }
    }
  }



template<typename eglue_type>
template<typename T1, typename T2>
arma_hot
inline
void
eglue_core<eglue_type>::apply_inplace_div(Cube<typename T1::elem_type>& out, const eGlueCube<T1, T2, eglue_type>& x)
  {
  arma_extra_debug_sigprint();
  
  const uword n_rows   = x.get_n_rows();
  const uword n_cols   = x.get_n_cols();
  const uword n_slices = x.get_n_slices();
  
  arma_debug_assert_same_size(out.n_rows, out.n_cols, out.n_slices, n_rows, n_cols, n_slices, "element-wise division");
  
  typedef typename T1::elem_type eT;
  
  eT* out_mem = out.memptr();
  
  const bool prefer_at_accessor = (ProxyCube<T1>::prefer_at_accessor || ProxyCube<T2>::prefer_at_accessor);
  
  if(prefer_at_accessor == false)
    {
    const bool  is_fixed = false;
    const uword n_elem   = out.n_elem;
    
    typename ProxyCube<T1>::ea_type P1 = x.P1.get_ea();
    typename ProxyCube<T2>::ea_type P2 = x.P2.get_ea();
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_1(/=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_1(/=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_1(/=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_1(/=, *); }
    }
  else
    {
    const ProxyCube<T1>& P1 = x.P1;
    const ProxyCube<T2>& P2 = x.P2;
    
         if(is_same_type<eglue_type, eglue_plus >::value == true) { arma_applier_3(/=, +); }
    else if(is_same_type<eglue_type, eglue_minus>::value == true) { arma_applier_3(/=, -); }
    else if(is_same_type<eglue_type, eglue_div  >::value == true) { arma_applier_3(/=, /); }
    else if(is_same_type<eglue_type, eglue_schur>::value == true) { arma_applier_3(/=, *); }
    }
  }



#undef arma_applier_1
#undef arma_applier_2
#undef arma_applier_3



//! @}
