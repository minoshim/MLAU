INCLUDE_PATH = $(HOME)/MLAU/include

# GNU compiler
CC = g++
CFLAGS= -O2 -march=native -fopenmp -lm -I$(INCLUDE_PATH)
# Intel compiler
# CC = icpc
# CFLAGS= -ipo -O2 -xhost -qopenmp -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd4c_1d.o prdc_1d.o

.PHONY : all
.PHONY : clean

all : $(OBJ)
	$(CC) -o a.out $(OBJ) $(CFLAGS)

%.o : %.c
	$(CC) -c $< $(CFLAGS)

%.o : %.cpp
	$(CC) -c $< $(CFLAGS)

clean : 
	rm *.o a.out dat/*.dat -f

mhd_fd4c_1d.o: prdc_1d.o $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h $(INCLUDE_PATH)/funcs.h
main.o: mhd_fd4c_1d.o prdc_1d.o myfunc.h $(INCLUDE_PATH)/funcs.h
