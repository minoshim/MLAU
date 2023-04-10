DIR0 = ../../..
include $(DIR0)/Makefile.inc
INCLUDE_PATH = $(DIR0)/include
INCDIR = -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd_ct_2d.o bc2d.o

.PHONY : all
.PHONY : clean

all : $(OBJ)
	$(CC) -o a.out $(OBJ) $(CFLAGS) $(INCDIR)

%.o : %.c
	$(CC) -c $< $(CFLAGS) $(INCDIR)

%.o : %.cpp
	$(CC) -c $< $(CFLAGS) $(INCDIR)

clean : 
	rm *.o a.out dat/*.dat -f

mhd_fd_ct_2d.o: bc2d.o $(INCLUDE_PATH)/funcs.h $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h
main.o: bc2d.o mhd_fd_ct_2d.o $(INCLUDE_PATH)/funcs.h myfunc.h cflcheck.hpp dataio.hpp global.hpp init.hpp new_delete.hpp setbc.hpp
