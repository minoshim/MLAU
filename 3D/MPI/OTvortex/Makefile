INCLUDE_PATH = $(HOME)/MLAU/include

CC = mpic++
# CFLAGS= -O2 -march=native -fopenmp -lm -I$(INCLUDE_PATH)
CFLAGS= -O2 -fopenmp -lm -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd_ct_3d.o mpi_sdrv3d.o mpi_bc3d.o

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

mhd_fd_ct_3d.o: mpi_sdrv3d.o mpi_bc3d.o $(INCLUDE_PATH)/funcs.h $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h

main.o: mhd_fd_ct_3d.o mpi_sdrv3d.o mpi_bc3d.o $(INCLUDE_PATH)/funcs.h myfunc.h bkup.hpp cflcheck.hpp dataio.hpp global.hpp init.hpp new_delete.hpp
