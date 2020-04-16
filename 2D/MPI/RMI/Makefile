INCLUDE_PATH = $(HOME)/MLAU/include

CC = mpic++
CFLAGS= -O2 -march=native -fopenmp -lm -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd_ct_2d.o mpi_sdrv2d.o mpi_sdrv2d_08.o mpi_xbc2d.o mpi_ybc2d.o

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

mhd_fd_ct_2d.o: mpi_sdrv2d.o mpi_sdrv2d_08.o mpi_xbc2d.o mpi_ybc2d.o $(INCLUDE_PATH)/funcs.h $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h

main.o: mhd_fd_ct_2d.o mpi_sdrv2d.o mpi_sdrv2d_08.o mpi_xbc2d.o mpi_ybc2d.o $(INCLUDE_PATH)/funcs.h myfunc.h bkup_load.hpp bkup_save.hpp cflcheck.hpp cflcomment.hpp dataio.hpp delete_vals.hpp global.hpp init_grid.hpp init_plasma.hpp new_vals.hpp inj_ro.hpp
