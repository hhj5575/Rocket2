# F90 Compiler for MPI
FC = mpif90

# Flags for Maximum Performance
FFLAGS = -O3
#FFLAGS = -O0 -debug full -check all -traceback -warn all

# Include Needed for Compile
IFLAGS = -I ../../include

# Library Needed for Linking
LFLAGS = -L ../../lib -ltecio -lstdc++

# OBJ List
OBJS = LIB_VTK_IO.o\
       converter.o

%.o: %.f90
	$(FC) -c $< -o $@ $(FFLAGS) $(IFLAGS)

all: converter

converter: $(OBJS)
	$(FC) -o $@ $(OBJS) $(FFLAGS) $(LFLAGS)

clean:
	find ./ -name '*~' -exec rm {} \;
	find ./ -name '*.o' -exec rm {} \;
	find ./ -name '*.mod' -exec rm {} \;
	find ./ -name '*__genmod.f90' -exec rm {} \;
	rm -f converter
