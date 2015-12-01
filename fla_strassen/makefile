#BLAS_LIB   := $(HOME)/flame/lib/libopenblas.a
FLAME_LIB  := $(HOME)/flame/lib/libflame.a
FLAME_INC  := $(HOME)/flame/include

# indicate where the object files are to be created
CC         := icc
FC         := ifort
LINKER     := $(CC)
#CFLAGS     := -O3 -Wall -I$(FLAME_INC) -msse3
#CFLAGS     := -O3 -Wall -I$(FLAME_INC) -g -pg
CFLAGS     := -O3 -Wall -I$(FLAME_INC)

SRC_PATH   := .
OBJ_PATH   := .

LDFLAGS       := -mkl=sequential -openmp -lm
#LDFLAGS       := -mkl=parallel -openmp -lm -g -pg

#LDFLAGS       := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -lcudart -lcublas -mkl=parallel -fopenmp -lgfortranbegin -lgfortran -lm
#LDFLAGS      := -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgfortran -lm -fopenmp
#LDFLAGS      := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lgfortran -lm -fopenmp -lcublas -lcudart
#LDFLAGS      := -I/usr/local/cuda/include -L/usr/local/cuda/lib64 -I/opt/intel/composer_xe_2013_sp1/mkl/include -L/opt/intel/composer_xe_2013_sp1/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -liomp5 -lmkl_core -lpthread -lgfortran -lm -fopenmp -lcublas -lcudart


# indicate where the FLAME include files reside
TEST_OBJS    := $(patsubst $(SRC_PATH)/%.c, $(OBJ_PATH)/%.o, $(wildcard $(SRC_PATH)/*.c))

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

driver.x: $(TEST_OBJS)
	$(LINKER) $(TEST_OBJS) $(FLAME_LIB) $(BLAS_LIB) $(LDFLAGS) -o driver.x

#test:   driver.x
#	echo "3 6.8 128 50 500 50" | ./driver.x 

clean:
	rm -f *.o *~ core *.x

