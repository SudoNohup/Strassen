#for devId in 1
#do
#  for size in $(seq 32 32 512)
#  do
#	#echo $devId
#	#echo $size
#	#nvcc -DSIZE=size dgemm.c -o dgemm -lcublas
#	#echo $devId $size
#	#./test_dgemm_cublas.x $devId $size
#	#./test_dsyrk_cublas.x $devId $size
#	./test_dtrsm_cublas.x $devId $size
#  done
#done

./driver.x 2000 500
./driver.x 2000 1000
./driver.x 2000 2000

#./driver.x 4000 500
#./driver.x 4000 1000
#./driver.x 4000 2000
#./driver.x 4000 4000

#./driver.x 8000 500
#./driver.x 8000 1000
#./driver.x 8000 2000
#./driver.x 8000 4000
#./driver.x 8000 8000
