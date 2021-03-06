#include "FLAME.h"
#include "Strassen_prototypes.h"

//
// --- FLA_Part_2x2() ----------------------------------------------------------
//
//

FLA_Error FLA_Part_Even_2x2_check( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                          FLA_Obj *A21, FLA_Obj *A22,
                              dim_t  mb,  dim_t     nb )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A22 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}



FLA_Error FLA_Part_Even_2x2( FLA_Obj *A,  FLA_Obj *A11, FLA_Obj *A12,
                                    FLA_Obj *A21, FLA_Obj *A22)
   //                     dim_t  mb,  dim_t     nb)
{
  FLA_Base_obj *base;
  dim_t         m, n, offm, offn;

//  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
//    FLA_Part_2x2_check( A,    A11, A12,
//                              A21, A22,    mb, nb);

  // Safeguard: if mb > m, reduce mb to m.
  //if ( mb > A.m ) mb = A.m;

  // Safeguard: if nb > n, reduce nb to n.
  //if ( nb > A.n ) nb = A.n;
  //

  m        = A->m;
  n        = A->n;
  offm     = A->offm;
  offn     = A->offn;
  base     = A->base;

  //printf("m:%d\t, n:%d\n", m, n);


  A11->m    = m / 2;
  A11->n    = n / 2;
  A11->offm = offm;
  A11->offn = offn;
  A11->base = base;

  A21->m    = m / 2;
  A21->n    = n / 2;
  A21->offm = offm + m / 2;
  A21->offn = offn;
  A21->base = base;

  A12->m    = m / 2;
  A12->n    = n / 2;
  A12->offm = offm;
  A12->offn = offn + n / 2;
  A12->base = base;

  A22->m    = m / 2;
  A22->n    = n / 2;
  A22->offm = offm + m / 2;
  A22->offn = offn + n / 2;
  A22->base = base;

  return FLA_SUCCESS;
}
