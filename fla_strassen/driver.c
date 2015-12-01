#include <stdio.h>
#include <math.h>
#include <time.h>

#include "FLAME.h"
#include "Strassen_prototypes.h"

int main(int argc, char **argv) {
  dim_t m, k, n, nb;
  int nrepeats, ireps;
  FLA_Obj A, B, C, Cref;
  double diff, dtime, dtime_best, gflops, dtime_average=0.0;

  FLA_Init();

  //n = 1024 * 4;
  //nb = 1024;

  if (argc != 5) {
	printf("argument number is wrong!\n");
	exit(0);
  }

  m  = atoi( argv[1] );
  k  = atoi( argv[2] );
  n  = atoi( argv[3] );

  nb = atoi( argv[4] );

  FLA_Obj_create( FLA_DOUBLE, m, k, 0, 0, &A );
  FLA_Obj_create( FLA_DOUBLE, k, n, 0, 0, &B );
  FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, &C );
  FLA_Obj_create( FLA_DOUBLE, m, n, 0, 0, &Cref );


  //FLA_Random_matrix( A );
  FLA_Set( FLA_ONE, A );
  //FLA_Random_matrix( B );
  FLA_Set( FLA_ONE, B );

  nrepeats = 3;
  gflops = 2.0 * m * k * n * 1.0e-9;

  /*
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
  */


  FLA_Set( FLA_ZERO, Cref );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, Cref );

  //for (nb = 128; nb <= 4096; nb <<= 1) {
  for (ireps = 0; ireps < nrepeats; ++ireps) {

	FLA_Set( FLA_ZERO, C);

	dtime = FLA_Clock();
	FLA_Strassen(A, B, C, nb);
	dtime = FLA_Clock() - dtime;

	if(ireps == 0)
	  dtime_best = dtime;
	else 
	  dtime_best = dtime < dtime_best ? dtime : dtime_best;

    dtime_average += dtime;
  }
  dtime_average /= nrepeats;
  //printf( "nb:%d\tFLA_Strassen:\ttime = [%le], flops = [%le]\n", nb, dtime_best, gflops / dtime_best);
  printf( "%lu\t%lu\t%le\t%le\t%le\t%le\n", n, nb, dtime_best, gflops / dtime_best, dtime_average, gflops / dtime_average);
  //}

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
