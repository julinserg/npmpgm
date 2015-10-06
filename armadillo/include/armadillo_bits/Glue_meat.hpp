// Copyright (C) 2008-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2008-2010 Conrad Sanderson
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


//! \addtogroup Glue
//! @{



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::Glue(const T1& in_A, const T2& in_B)
  : A(in_A)
  , B(in_B)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::Glue(const T1& in_A, const T2& in_B, const uword in_aux_uword)
  : A(in_A)
  , B(in_B)
  , aux_uword(in_aux_uword)
  {
  arma_extra_debug_sigprint();
  }



template<typename T1, typename T2, typename glue_type>
inline
Glue<T1,T2,glue_type>::~Glue()
  {
  arma_extra_debug_sigprint();
  }



//! @}
