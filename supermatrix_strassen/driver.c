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

  FLASH_Queue_set_num_threads( 24 );
  FLASH_Queue_set_sorting( 0 );
  FLASH_Queue_set_caching( 0 );
  FLASH_Queue_set_work_stealing( 0 );
  FLASH_Queue_set_data_affinity( 0 );
  //FLASH_Queue_set_verbose_output( FLASH_QUEUE_VERBOSE_NONE);
  //FLASH_Queue_set_verbose_output( FLASH_QUEUE_VERBOSE_READABLE);
  //FLASH_Queue_set_verbose_output( FLASH_QUEUE_VERBOSE_ELSE);

  n = 1024 * 4;
  //nb_alg = 1024;
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &A );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &B );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &C );
  FLA_Obj_create( FLA_DOUBLE, n, n, 0, 0, &Cref );

  FLA_Random_matrix( A );
  //FLA_Set( FLA_ONE, A );
  FLA_Random_matrix( B );
  //FLA_Set( FLA_ONE, B );

  nrepeats = 5;
  gflops = 2.0 * n * n * n * 1.0e-9;

  for (ireps = 0; ireps < nrepeats; ++ireps) {

	FLA_Set( FLA_ZERO, Cref );

	dtime = FLA_Clock();
	printf("start FLA_Gemm\n");
	FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ONE, Cref );
	printf("finish FLA_Gemm\n");
	fflush(stdout);
	dtime = FLA_Clock() - dtime;

	if(ireps == 0)
	  dtime_best = dtime;
	else 
	  dtime_best = dtime < dtime_best ? dtime : dtime_best;

	//FLA_Strassen(A, B, C);
  }
  printf( "FLA_Gemm:\ttime = [%le], flops = [%le]\n", dtime_best, gflops / dtime_best);


  //FLA_Obj_show("Cref:", Cref, "%11.3e", "...END...");
  //FLA_Obj_show("C:", C, "%11.3e", "...END...");

  //for (nb_alg = 1024; nb_alg <= 1024; nb_alg <<= 1) {
  for (nb_alg = 128; nb_alg <= 4096; nb_alg <<= 1) {
  for (ireps = 0; ireps < nrepeats; ++ireps) {

	FLA_Set( FLA_ZERO, C);

	dtime = FLA_Clock();

	FLASH_Queue_begin();

	Strassen_Workspace *wks = Strassen_Workspace_new();

	FLASH_Obj_create_without_buffer( FLA_Obj_datatype(A), FLA_Obj_length(A), FLA_Obj_width(A), 1, &nb_alg, &AH );
	FLASH_Obj_attach_buffer(FLA_Obj_buffer_at_view(A), FLA_Obj_row_stride(A), FLA_Obj_col_stride(A), &AH);
	FLASH_Obj_create_without_buffer( FLA_Obj_datatype(B), FLA_Obj_length(B), FLA_Obj_width(B), 1, &nb_alg, &BH );
	FLASH_Obj_attach_buffer(FLA_Obj_buffer_at_view(B), FLA_Obj_row_stride(B), FLA_Obj_col_stride(B), &BH);
	FLASH_Obj_create_without_buffer( FLA_Obj_datatype(C), FLA_Obj_length(C), FLA_Obj_width(C), 1, &nb_alg, &CH );
	FLASH_Obj_attach_buffer(FLA_Obj_buffer_at_view(C), FLA_Obj_row_stride(C), FLA_Obj_col_stride(C), &CH);

	FLASH_Strassen(&AH, &BH, &CH, wks, nb_alg);

	
	//Strassen_Workspace_print(wks);

	//printf("before execution stage!\n");

	FLASH_Queue_end();

	//printf("out of execution stage!\n");
	dtime = FLA_Clock() - dtime;
	Strassen_Workspace_free(wks);

	//printf("free workspace!\n");

	FLASH_Obj_free_without_buffer( &AH );
	FLASH_Obj_free_without_buffer( &BH );
	FLASH_Obj_free_without_buffer( &CH );


	if(ireps == 0)
	  dtime_best = dtime;
	else 
	  dtime_best = dtime < dtime_best ? dtime : dtime_best;
  }
  diff = FLA_Max_elemwise_diff( C, Cref );
  //printf( "diff = [ %le]\n", diff );

  printf( "nb_alg:%lu\tFLA_Strassen:\ttime = [%le], flops = [%le], diff = [%le]\n", nb_alg, dtime_best, gflops / dtime_best, diff);
  }


  //FLA_Obj_show("C:", C, "%11.3e", "...END...");
  //FLA_Obj_show("Cref:", Cref, "%11.3e", "...END...");
  //diff = FLA_Max_elemwise_diff( C, Cref );
  //printf( "diff = [ %le]\n", diff );
  fflush( stdout );


  FLA_Obj_free( &A );
  FLA_Obj_free( &B );
  FLA_Obj_free( &C );
  FLA_Obj_free( &Cref );

  FLA_Finalize();

  return 0;
}
