#BLAS_LIB   := $(HOME)/flame/lib/libopenblas.a
FLAME_LIB  := $(HOME)/flame/lib/libflame.a
FLAME_INC  := $(HOME)/flame/include

# indicate where the object files are to be created
#CC         := gcc
CC         := icc
FC         := gfortran
LINKER     := $(CC)
#CFLAGS     := -O3 -Wall -I$(FLAME_INC) -msse3
CFLAGS     := -O3 -I$(FLAME_INC) -msse3
FFLAGS     := $(CFLAGS)

SRC_PATH   := .
OBJ_PATH   := .

#LDFLAGS      := -L/usr/lib/gcc/x86_64-linux-gnu/4.7 -L/usr/lib/gcc/x86_64-linux-gnu/4.7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.7/../../.. -lgfortranbegin -lgfortran -lquadmath -lpthread -lm
#LDFLAGS       := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcublas -mkl=sequential -fopenmp -lgfortranbegin -lgfortran -lm
#LDFLAGS       := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcublas -mkl=parallel -fopenmp -lgfortranbegin -lgfortran -lm
#LDFLAGS      := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgfortran -lm -fopenmp -lcublas -lcudart
#LDFLAGS      := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lpthread -lgfortran -lm -fopenmp -lcublas -lcudart
#LDFLAGS      := -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lgfortran -lm -fopenmp

LDFLAGS      := -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lgfortran -lm -fopenmp
#LDFLAGS      := -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lgfortran -lm -pthread
#LDFLAGS      := -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lgfortran -lm -fopenmp

# indicate where the FLAME include files reside

#TEST_OBJS  := driver.o \
#	Gemm_unb_var1.o Gemm_blk_var1.o \
#	Gemm_unb_var2.o Gemm_blk_var2.o \
#	Gemm_unb_var3.o Gemm_blk_var3.o \
#	Gepp_blk_var1.o Gebp.o

TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

%.o: %.c
	echo $(TEST_OBJS)
	$(CC) $(CFLAGS) -c $< -o $@
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)
	#$(LINKER) $(TEST_OBJS) $(FLAME_LIB) $(BLAS_LIB) -llapack $(LDFLAGS) -o driver.x
	$(LINKER) $(TEST_OBJS) $(FLAME_LIB) $(BLAS_LIB) $(LDFLAGS) -o driver.x

test:   driver.x
	echo "3 6.8 128 50 500 50" | ./driver.x 
clean:
	rm -f *.o *~ core *.x

