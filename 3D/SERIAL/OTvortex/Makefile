INCLUDE_PATH = $(HOME)/MLAU/include

# GNU compiler
CC = g++
CFLAGS= -O2 -march=native -fopenmp -lm -I$(INCLUDE_PATH)
# Intel compiler
# CC = icpc
# CFLAGS= -ipo -O2 -xhost -qopenmp -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd_ct_3d.o bc3d.o

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

mhd_fd_ct_3d.o: bc3d.o $(INCLUDE_PATH)/funcs.h $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h

main.o: bc3d.o mhd_fd_ct_3d.o $(INCLUDE_PATH)/funcs.h myfunc.h cflcheck.hpp cflcomment.hpp dataio.hpp delete_vals.hpp global.hpp init_grid.hpp init_plasma.hpp new_vals.hpp setbc.hpp
