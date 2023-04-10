DIR0 = ../../..
include $(DIR0)/Makefile.inc
INCLUDE_PATH = $(DIR0)/include
INCDIR = -I$(INCLUDE_PATH)

OBJ = main.o mhd_fd4c_1d.o prdc_1d.o

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

mhd_fd4c_1d.o: prdc_1d.o $(INCLUDE_PATH)/interp.h $(INCLUDE_PATH)/mhd_flux.h $(INCLUDE_PATH)/mhd_eigen.h $(INCLUDE_PATH)/funcs.h
main.o: mhd_fd4c_1d.o prdc_1d.o myfunc.h $(INCLUDE_PATH)/funcs.h
