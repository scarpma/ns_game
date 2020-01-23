FFTW3 = /home/scarpma/fftw-3.3.8/.libs/libfftw3.a
FFTW3PATH = /home/scarpma/fftw-3.3.8/api
FC = gfortran#~/intel/compilers_and_libraries_2020/linux/bin/intel64/ifort
FLAGS = -O3 -fmax-errors=3 -ffree-line-length-512 #-fdec-structure -fbounds-check
TARGET = ns

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)
PATHOBJS = $(OBJS:%=objs/%)

objs/%.o: %.f90 
	$(FC) $(FLAGS) -c $< -o $@ -Jobjs -I$(FFTW3PATH)

$(TARGET): $(PATHOBJS)
	$(FC) -o $(TARGET) $(PATHOBJS) $(LIBS) $(FFTW3)

clean:
	rm objs/*.o objs/*.mod $(TARGET)

include make.dep
