
FC = gfortran
FLAGS = -O2
TARGET = sfs

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm *.o *.mod

distclean:
	rm $(TARGET)


main.o : precision.o  common.o inout.o
common.o : linear_solver.o precision.o mat_def.o
inout.o: precision.o
linear_solver.o: precision.o mat_def.o
mat_def.o: allocation.o precision.o
allocation.o: precision.o
