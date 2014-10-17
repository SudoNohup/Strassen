#include "FLAME.h"
#include "Strassen_prototypes.h"

//FLA_Error FLA_Part_even_2x2( FLA_Obj A, FLA_Obj *A11, FLA_Obj *A12, FLA_Obj *A21, FLA_Obj *A22, dim_t mb, dim_t nb
//


extern fla_copy_t* flash_copy_cntl;
extern fla_axpy_t* flash_axpy_cntl;
extern fla_gemm_t* flash_gemm_cntl_mm_op;

//Workspace...
FLA_Error FLASH_Strassen(FLA_Obj AH, FLA_Obj BH, FLA_Obj CH, Strassen_Workspace *wks, dim_t nb_alg) {
  int n;
  FLA_Obj A11H, A12H, A21H, A22H;
  FLA_Obj B11H, B12H, B21H, B22H;
  FLA_Obj C11H, C12H, C21H, C22H;
  FLA_Obj S1, S2, S3, S4, S5, S6, S7;
  FLA_Obj S1H, S2H, S3H, S4H, S5H, S6H, S7H;
  FLA_Obj T1, T2, T3, T4, T5, T6, T7;
  FLA_Obj T1H, T2H, T3H, T4H, T5H, T6H, T7H;
  FLA_Obj M1, M2, M3, M4, M5, M6, M7;
  FLA_Obj M1H, M2H, M3H, M4H, M5H, M6H, M7H;

  //we need to check the dimension of A, B, C ?? (mxn) * (nxk) = (m*k);

  // Resursive base
  // A's dimension is less than 128, assume this is a square matrix....
  /*
  if ( FLA_Obj_length(A) <= nb_alg ) {
	FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, C);
	//printf("C(%d %d)", C.offm, C.offn);
	//FLA_Obj_show("A:", A, "%11.3e", "...END...");
	//FLA_Obj_show("B:", B, "%11.3e", "...END...");
	//FLA_Obj_show("C:", C, "%11.3e", "...END...");
	return FLA_SUCCESS;
  }
  */

  if ( FLA_Obj_length(AH) <= 1 ) {
	printf("Enter the base of recursion!\n");
	FLA_Gemm_internal( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, AH, BH, FLA_ONE, CH, flash_gemm_cntl_mm_op);
	//printf("C(%d %d)", C.offm, C.offn);
	//FLA_Obj_show("A:", A, "%11.3e", "...END...");
	//FLA_Obj_show("B:", B, "%11.3e", "...END...");
	//FLA_Obj_show("C:", C, "%11.3e", "...END...");
	return FLA_SUCCESS;
  }


  //n = FLA_Obj_length( *FLASH_OBJ_PTR_AT( A ) ) / 2;
  n = FLA_Obj_length( AH ) / 2 * nb_alg;


  printf("S1 addr: %d\n", &S1);
  printf("S1 addr: %d\n", &S2);
  printf("S1H addr: %d\n", &S1H);

  Strassen_allocateObj(wks, &S1, &S1H, n, nb_alg);
  Strassen_allocateObj(wks, &S2, &S2H, n, nb_alg);
  Strassen_allocateObj(wks, &S3, &S3H, n, nb_alg);
  Strassen_allocateObj(wks, &S4, &S4H, n, nb_alg);
  Strassen_allocateObj(wks, &S5, &S5H, n, nb_alg);
  Strassen_allocateObj(wks, &S6, &S6H, n, nb_alg);
  Strassen_allocateObj(wks, &S7, &S7H, n, nb_alg);
  Strassen_allocateObj(wks, &T1, &T1H, n, nb_alg);
  Strassen_allocateObj(wks, &T2, &T2H, n, nb_alg);
  Strassen_allocateObj(wks, &T3, &T3H, n, nb_alg);
  Strassen_allocateObj(wks, &T4, &T4H, n, nb_alg);
  Strassen_allocateObj(wks, &T5, &T5H, n, nb_alg);
  Strassen_allocateObj(wks, &T6, &T6H, n, nb_alg);
  Strassen_allocateObj(wks, &T7, &T7H, n, nb_alg);
  Strassen_allocateObj(wks, &M1, &M1H, n, nb_alg);
  Strassen_allocateObj(wks, &M2, &M2H, n, nb_alg);
  Strassen_allocateObj(wks, &M3, &M3H, n, nb_alg);
  Strassen_allocateObj(wks, &M4, &M4H, n, nb_alg);
  Strassen_allocateObj(wks, &M5, &M5H, n, nb_alg);
  Strassen_allocateObj(wks, &M6, &M6H, n, nb_alg);
  Strassen_allocateObj(wks, &M7, &M7H, n, nb_alg);

  FLA_Part_Even_2x2( AH,  &A11H, &A12H,
                         &A21H, &A22H);

  FLA_Part_Even_2x2( BH,  &B11H, &B12H,
                         &B21H, &B22H);

  FLA_Part_Even_2x2( CH,  &C11H, &C12H,
	                     &C21H, &C22H);


  printf("A11H_m:%d\n", FLA_Obj_length(A11H));
  printf("S1H_m:%d\n", FLA_Obj_length(S1H));
  printf("A11H_n:%d\n:", FLA_Obj_width(A11H));
  printf("S1H_n:%d\n", FLA_Obj_width(S1H));

  printf("flag1\n");
  FLA_Copy_internal(A11H, S1H, flash_copy_cntl);

  FLASH_print_struct(A11H);
  FLASH_print_struct(S1H);

  FLA_Obj_show("A11:", *FLASH_OBJ_PTR_AT(A11H), "%11.3e", "...END...");
  FLA_Obj_show("S1:", S1, "%11.3e", "...END...");
  //We need to verify the Copy result here...

  printf("flag2\n");
  FLA_Axpy_internal(FLA_ONE, A22H, S1H, flash_axpy_cntl);
  printf("flag3\n");


  FLA_Copy_internal(A21H, S2H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, A22H, S2H, flash_axpy_cntl);

  FLA_Copy_internal(A11H, S3H, flash_copy_cntl);

  FLA_Copy_internal(A22H, S4H, flash_copy_cntl);

  FLA_Copy_internal(A11H, S5H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, A12H, S5H, flash_axpy_cntl);

  FLA_Copy_internal(A21H, S6H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, A11H, S6H, flash_axpy_cntl);

  FLA_Copy_internal(A12H, S7H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, A22H, S7H, flash_axpy_cntl);

  FLA_Copy_internal(B11H, T1H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, B22H, T1H, flash_axpy_cntl);

  FLA_Copy_internal(B11H, T2H, flash_copy_cntl);

  FLA_Copy_internal(B12H, T3H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, B22H, T3H, flash_axpy_cntl);

  FLA_Copy_internal(B21H, T4H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, B11H, T4H, flash_axpy_cntl);

  FLA_Copy_internal(B22H, T5H, flash_copy_cntl);

  FLA_Copy_internal(B11H, T6H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, B12H, T6H, flash_axpy_cntl);

  FLA_Copy_internal(B21H, T7H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, B22H, T7H, flash_axpy_cntl);


  printf("flag5\n");

  FLASH_Strassen(S1H, T1H, M1H, wks, nb_alg);
  FLASH_Strassen(S2H, T2H, M2H, wks, nb_alg);
  FLASH_Strassen(S3H, T3H, M3H, wks, nb_alg);
  FLASH_Strassen(S4H, T4H, M4H, wks, nb_alg);
  FLASH_Strassen(S5H, T5H, M5H, wks, nb_alg);
  FLASH_Strassen(S6H, T6H, M6H, wks, nb_alg);
  FLASH_Strassen(S7H, T7H, M7H, wks, nb_alg);

  printf("flag6\n");

  FLA_Copy_internal(M1H, C11H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, M4H, C11H, flash_axpy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, M5H, C11H, flash_axpy_cntl);
  FLA_Axpy_internal(FLA_ONE, M7H, C11H, flash_axpy_cntl);

  printf("flag7\n");
  FLA_Copy_internal(M3H, C12H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, M5H, C12H, flash_axpy_cntl);

  FLA_Copy_internal(M2H, C21H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_ONE, M4H, C21H, flash_axpy_cntl);

  printf("flag8\n");
  FLA_Copy_internal(M1H, C22H, flash_copy_cntl);
  FLA_Axpy_internal(FLA_MINUS_ONE, M2H, C22H, flash_axpy_cntl);
  FLA_Axpy_internal(FLA_ONE, M3H, C22H, flash_axpy_cntl);
  FLA_Axpy_internal(FLA_ONE, M6H, C22H, flash_axpy_cntl);

  return FLA_SUCCESS;

}



