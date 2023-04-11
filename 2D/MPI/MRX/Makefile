DIR0 = ../../..
include $(DIR0)/Makefile.inc
INCLUDE_PATH = $(DIR0)/include
INCDIR = -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd_ct_2d.o mpi_sdrv2d.o mpi_bc2d.o

.PHONY : all
.PHONY : clean

all : $(OBJ)
	$(MPICC) -o a.out $(OBJ) $(CFLAGS) $(INCDIR)

%.o : %.c
	$(MPICC) -c $< $(CFLAGS) $(INCDIR)

%.o : %.cpp
	$(MPICC) -c $< $(CFLAGS) $(INCDIR)

clean : 
	rm *.o a.out dat/*.dat -f

mhd_fd_ct_2d.o: mpi_sdrv2d.o mpi_bc2d.o $(INCLUDE_PATH)/funcs.h $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h
main.o: mhd_fd_ct_2d.o mpi_sdrv2d.o mpi_bc2d.o $(INCLUDE_PATH)/funcs.h myfunc.h bkup.hpp cflcheck.hpp dataio.hpp global.hpp init.hpp new_delete.hpp
