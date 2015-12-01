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

#./driver.x 2000 500
#./driver.x 2000 1000
#./driver.x 2000 2000

#./driver.x 4000 500
#./driver.x 4000 1000
#./driver.x 4000 2000
#./driver.x 4000 4000

#./driver.x 8000 500
#./driver.x 8000 1000
#./driver.x 8000 2000
#./driver.x 8000 4000
#./driver.x 8000 8000


#./driver.x 4000 2000
#./driver.x 8000 4000

##for (( size=128; size<1025; size+=128 ))
#for (( size=512; size<4500; size+=512 ))
#do
#  size2=$((size/2))
#  echo $size2
#  ./driver.x $size $size2
#done


#./driver.x 4096 512
#./driver.x 4096 1024
#./driver.x 4096 2048
#./driver.x 4096 4096



#./driver.x 128 128 128 64
#./driver.x 128 32 128 64
#./driver.x 4096 4096 4096 2048
#./driver.x 4096 512 4096 2048
#./driver.x 4096 512 4096 4096
#./driver.x 4096 4096 4096 4096
#./driver.x 4096 256 4096 4096
#./driver.x 4096 32 4096 4096
./driver.x 4096 16 4096 4096






