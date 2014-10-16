#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"
#include "Strassen_prototypes.h"


int main() {
  dim_t n, nb_alg;
  int nrepeats, ireps;
  FLA_Obj A, B, C, Cref;
  FLA_Obj AH, BH, CH;
  double diff, dtime, dtime_best, gflops;

  FLA_Init();

  n = 128;
  nb_alg = 4;
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &A );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &B );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &C );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &Cref );

  //FLA_Random_matrix( A );
  FLA_Set( FLA_ONE, A );
  //FLA_Random_matrix( B );
  FLA_Set( FLA_ONE, B );

  nrepeats = 1;
  gflops = 2.0 * n * n * n * 1.0e-9;

  for (ireps = 0; ireps < nrepeats; ++ireps) {

	FLA_Set( FLA_ZERO, Cref );

	dtime = FLA_Clock();
	FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, Cref );
	dtime = FLA_Clock() - dtime;

	if(ireps == 0)
	  dtime_best = dtime;
	else 
	  dtime_best = dtime < dtime_best ? dtime : dtime_best;

	//FLA_Strassen(A, B, C);
  }

  printf( "FLA_Gemm:\ttime = [%le], flops = [%le]\n", dtime_best, gflops / dtime_best);


  for (ireps = 0; ireps < nrepeats; ++ireps) {

	FLA_Set( FLA_ZERO, C);

	dtime = FLA_Clock();

	FLASH_Queue_begin();

	Strassen_Workspace *wks = Strassen_Workspace_new();

	FLASH_Obj_create_without_buffer( FLA_Obj_datetype(A), FLA_Obj_length(A), FLA_Obj_width(A), 1, &nb_alg, &AH );
	FLASH_attach_buffer(FLA_Obj_buffer_at_view(A), FLA_Obj_row_stride(A), FLA_Obj_col_stride(A), &AH);
	FLASH_Obj_create_without_buffer( FLA_Obj_datetype(B), FLA_Obj_length(B), FLA_Obj_width(B), 1, &nb_alg, &BH );
	FLASH_attach_buffer(FLA_Obj_buffer_at_view(B), FLA_Obj_row_stride(B), FLA_Obj_col_stride(B), &BH);
	FLASH_Obj_create_without_buffer( FLA_Obj_datetype(C), FLA_Obj_length(C), FLA_Obj_width(C), 1, &nb_alg, &CH );
	FLASH_attach_buffer(FLA_Obj_buffer_at_view(C), FLA_Obj_row_stride(C), FLA_Obj_col_stride(C), &CH);

	FLASH_Strassen(A, B, C, wks, nb_alg);
	FLASH_Queue_end();
	dtime = FLA_Clock() - dtime;

	if(ireps == 0)
	  dtime_best = dtime;
	else 
	  dtime_best = dtime < dtime_best ? dtime : dtime_best;
  }


  printf( "FLA_Strassen:\ttime = [%le], flops = [%le]\n", dtime_best, gflops / dtime_best);

  diff = FLA_Max_elemwise_diff( C, Cref );
  printf( "diff = [ %le]\n", diff );
  fflush( stdout );

  FLA_Obj_free( &A );
  FLA_Obj_free( &B );
  FLA_Obj_free( &C );
  FLA_Obj_free( &Cref );

  FLA_Finalize();

  return 0;
}