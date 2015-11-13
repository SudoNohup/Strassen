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

#./driver.x 4096 512 16 1
#./driver.x 4096 512 8 2
#./driver.x 4096 512 4 4
#./driver.x 4096 512 2 8
#./driver.x 4096 512 1 16
#
#./driver.x 4096 1024 16 1
#./driver.x 4096 1024 8 2
#./driver.x 4096 1024 4 4
#./driver.x 4096 1024 2 8
#./driver.x 4096 1024 1 16
#
#./driver.x 4096 2048 16 1
#./driver.x 4096 2048 8 2
#./driver.x 4096 2048 4 4
#./driver.x 4096 2048 2 8
#./driver.x 4096 2048 1 16

#./driver.x 4096 4096 16 1
#./driver.x 4096 4096 8 2
#./driver.x 4096 4096 4 4
#./driver.x 4096 4096 2 8
#./driver.x 4096 4096 1 16

#./driver.x 2000 500 16 1
#./driver.x 2000 500 8 2
#./driver.x 2000 500 4 4
#./driver.x 2000 500 2 8
#./driver.x 2000 500 1 16

./driver.x 2000 1000 16 1
./driver.x 2000 1000 8 2
./driver.x 2000 1000 4 4
./driver.x 2000 1000 2 8
./driver.x 2000 1000 1 16

