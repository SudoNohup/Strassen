#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"
#include "Strassen_prototypes.h"


int main() {
  dim_t n;
  int nrepeats, ireps;
  FLA_Obj A, B, C, Cref;
  double diff, dtime, dtime_best, gflops;

  FLA_Init();

  n = 1024 * 4;
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
	FLA_Strassen(A, B, C);
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
