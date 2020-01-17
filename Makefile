FC = gfortran#~/intel/compilers_and_libraries_2020/linux/bin/intel64/ifort
FLAGS = -O3 -fmax-errors=3 -ffree-line-length-512 #-fdec-structure -fbounds-check
TARGET = ns

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)
PATHOBJS = $(OBJS:%=objs/%)

objs/%.o: %.f90 
	$(FC) $(FLAGS) -c $< -o $@ -Jobjs

$(TARGET): $(PATHOBJS)
	$(FC) -o $(TARGET) $(PATHOBJS) $(LIBS)

clean:
	rm objs/*.o objs/*.mod $(TARGET)

include make.dep
