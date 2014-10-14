#include "FLAME.h"
#include "Strassen_prototypes.h"

//FLA_Error FLA_Part_even_2x2( FLA_Obj A, FLA_Obj *A11, FLA_Obj *A12, FLA_Obj *A21, FLA_Obj *A22, dim_t mb, dim_t nb


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


//
// --- FLA_Part_2x2() ----------------------------------------------------------
//

FLA_Error FLA_Part_Even_2x2( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
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

  m        = A.m;
  n        = A.n;
  offm     = A.offm;
  offn     = A.offn;
  base     = A.base;


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


FLA_Error FLA_Strassen(FLA_Obj A, FLA_Obj B, FLA_Obj C) {
  int n;
  FLA_Obj A11, A12, A21, A22;
  FLA_Obj B11, B12, B21, B22;
  FLA_Obj C11, C12, C21, C22;
  FLA_Obj S1, S2, S3, S4, S5, S6, S7;
  FLA_Obj T1, T2, T3, T4, T5, T6, T7;
  FLA_Obj M1, M2, M3, M4, M5, M6, M7;

  //we need to check the dimension of A, B, C ?? (mxn) * (nxk) = (m*k);

  // Resursive base
  // A's dimension is less than 128, assume this is a square matrix....
  if ( FLA_Obj_length(A) <= 256 ) {
	FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, C);
	//printf("C(%d %d)", C.offm, C.offn);
	//FLA_Obj_show("A:", A, "%11.3e", "...END...");
	//FLA_Obj_show("B:", B, "%11.3e", "...END...");
	//FLA_Obj_show("C:", C, "%11.3e", "...END...");
	return FLA_SUCCESS;
  }

  n = FLA_Obj_length( A ) / 2;
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S1 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S2 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S3 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S4 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S5 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S6 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &S7 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T1 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T2 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T3 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T4 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T5 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T6 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &T7 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M1 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M2 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M3 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M4 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M5 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M6 );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &M7 );


  FLA_Set( FLA_ZERO, S1 );
  FLA_Set( FLA_ZERO, S2 );
  FLA_Set( FLA_ZERO, S3 );
  FLA_Set( FLA_ZERO, S4 );
  FLA_Set( FLA_ZERO, S5 );
  FLA_Set( FLA_ZERO, S6 );
  FLA_Set( FLA_ZERO, S7 );
  FLA_Set( FLA_ZERO, T1 );
  FLA_Set( FLA_ZERO, T2 );
  FLA_Set( FLA_ZERO, T3 );
  FLA_Set( FLA_ZERO, T4 );
  FLA_Set( FLA_ZERO, T5 );
  FLA_Set( FLA_ZERO, T6 );
  FLA_Set( FLA_ZERO, T7 );
  FLA_Set( FLA_ZERO, M1 );
  FLA_Set( FLA_ZERO, M2 );
  FLA_Set( FLA_ZERO, M3 );
  FLA_Set( FLA_ZERO, M4 );
  FLA_Set( FLA_ZERO, M5 );
  FLA_Set( FLA_ZERO, M6 );
  FLA_Set( FLA_ZERO, M7 );



  FLA_Part_Even_2x2( A,  &A11, &A12,
                         &A21, &A22);

  FLA_Part_Even_2x2( B,  &B11, &B12,
                         &B21, &B22);

  FLA_Part_Even_2x2( C,  &C11, &C12,
	                     &C21, &C22);

  FLA_Copy(A11, S1);
  FLA_Axpy(FLA_ONE, A22, S1); 
  //FLA_Obj_show("S1:", S1, "%11.3e", "...END...");

  FLA_Copy(A21, S2);
  FLA_Axpy(FLA_ONE, A22, S2); 

  FLA_Copy(A11, S3);

  FLA_Copy(A22, S4);

  FLA_Copy(A11, S5);
  FLA_Axpy(FLA_ONE, A12, S5); 

  FLA_Copy(A21, S6);
  FLA_Axpy(FLA_MINUS_ONE, A11, S6); 

  FLA_Copy(A12, S7);
  FLA_Axpy(FLA_MINUS_ONE, A22, S7); 

  
  FLA_Copy(B11, T1);
  FLA_Axpy(FLA_ONE, B22, T1);

  FLA_Copy(B11, T2);

  FLA_Copy(B12, T3);
  FLA_Axpy(FLA_MINUS_ONE, B22, T3);

  FLA_Copy(B21, T4);
  FLA_Axpy(FLA_MINUS_ONE, B11, T4);

  FLA_Copy(B22, T5);
 
  FLA_Copy(B11, T6);
  FLA_Axpy(FLA_ONE, B12, T6);

  FLA_Copy(B21, T7);
  FLA_Axpy(FLA_ONE, B22, T7);

  FLA_Strassen(S1, T1, M1);

  //FLA_Obj_show("S1:", S1, "%11.3e", "...END...");
  //FLA_Obj_show("T1:", T1, "%11.3e", "...END...");
  //FLA_Obj_show("M1:", M1, "%11.3e", "...END...");
  //printf("----------------------------------------------------S1, T1-------------------M1\n");
  FLA_Strassen(S2, T2, M2);
  //FLA_Obj_show("S2:", S2, "%11.3e", "...END...");
  //FLA_Obj_show("T2:", T2, "%11.3e", "...END...");
  //FLA_Obj_show("M2:", M2, "%11.3e", "...END...");
  //printf("----------------------------------------------------S2, T2-------------------M2\n");
  //FLA_Obj_show("M3:", M3, "%11.3e", "...END...");
 
  FLA_Strassen(S3, T3, M3);
  FLA_Strassen(S4, T4, M4);
  FLA_Strassen(S5, T5, M5);
  //FLA_Obj_show("S5:", S5, "%11.3e", "...END...");
  //FLA_Obj_show("T5:", T5, "%11.3e", "...END...");
  //FLA_Obj_show("M5:", M5, "%11.3e", "...END...");
  //printf("----------------------------------------------------S5, T5-------------------M2\n");
 
  FLA_Strassen(S6, T6, M6);
  FLA_Strassen(S7, T7, M7);

  FLA_Copy(M1, C11);
  FLA_Axpy(FLA_ONE, M4, C11);
  FLA_Axpy(FLA_MINUS_ONE, M5, C11);
  FLA_Axpy(FLA_ONE, M7, C11);

  FLA_Copy(M3, C12);
  FLA_Axpy(FLA_ONE, M5, C12);

  FLA_Copy(M2, C21);
  FLA_Axpy(FLA_ONE, M4, C21);

  FLA_Copy(M1, C22);
  FLA_Axpy(FLA_MINUS_ONE, M2, C22);
  FLA_Axpy(FLA_ONE, M3, C22);
  FLA_Axpy(FLA_ONE, M6, C22);

  FLA_Obj_free( &S1 );
  FLA_Obj_free( &S2 );
  FLA_Obj_free( &S3 );
  FLA_Obj_free( &S4 );
  FLA_Obj_free( &S5 );
  FLA_Obj_free( &S6 );
  FLA_Obj_free( &S7 );
  FLA_Obj_free( &T1 );
  FLA_Obj_free( &T2 );
  FLA_Obj_free( &T3 );
  FLA_Obj_free( &T4 );
  FLA_Obj_free( &T5 );
  FLA_Obj_free( &T6 );
  FLA_Obj_free( &T7 );
  FLA_Obj_free( &M1 );
  FLA_Obj_free( &M2 );
  FLA_Obj_free( &M3 );
  FLA_Obj_free( &M4 );
  FLA_Obj_free( &M5 );
  FLA_Obj_free( &M6 );
  FLA_Obj_free( &M7 );


  return FLA_SUCCESS;

}




